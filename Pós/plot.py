import numpy as np
import pyvista as pv
import re


def ler_blocos_resultados(filepath):
    """
    Lê todos os blocos de resultados do arquivo de texto (resultados FEM),
    incluindo Tensões, Deformações nos Centros e nos Pontos de Gauss (GP1 a GP9),
    além de deslocamentos nodais.

    Retorna:
        nodes: {id: [x, y, 0.0]}
        elements: {id: [n1, n2, ...]} (0-based)
        displacements: {id: [dx, dy, 0.0]}
        dados_por_bloco: dict do tipo:
            {'Tensoes_Centros dos Elementos': {id_elem: [xx, yy, xy]}, 
             'Deformacoes_GP 2': {...}, etc.}
    """
    nodes = {}
    elements = {}
    displacements = {}
    dados_por_bloco = {}

    current_section = None
    current_bloco_nome = None

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Detecta cabeçalhos
            if line.startswith("# Coordenadas dos Nos"):
                current_section = 'nodes'
                current_bloco_nome = None # Reset bloco nome
                continue
            elif line.startswith("# Conectividade dos Elementos"):
                current_section = 'elements'
                current_bloco_nome = None # Reset bloco nome
                continue
            elif line.startswith("# Deslocamentos Nodais"):
                current_section = 'displacements'
                current_bloco_nome = None # Reset bloco nome
                continue
            elif line.startswith("# Tensoes") or line.startswith("# Deformacoes"):
                match = re.match(r"# (Tensoes|Deformacoes) (.*)", line)
                if match:
                    tipo, local_full = match.groups()
                    # Extract only the relevant part of 'local_full'
                    if "nos Centros dos Elementos" in local_full:
                        local = "Centros dos Elementos"
                    elif "no GP" in local_full:
                        local_match = re.search(r"no (GP \d+)", local_full)
                        if local_match:
                            local = local_match.groups()[0]
                        else:
                            local = local_full.split('(')[0].strip() # Fallback
                    else:
                        local = local_full.split('(')[0].strip() # Fallback for other cases

                    current_section = 'campo'
                    current_bloco_nome = f"{tipo}_{local}"
                    dados_por_bloco[current_bloco_nome] = {}
                continue
            elif line.startswith("#"):
                continue  # outros comentários

            parts = line.split()
            if current_section == 'nodes' and len(parts) >= 3:
                node_id = int(parts[0])
                nodes[node_id] = [float(parts[1]), float(parts[2]), 0.0]
            elif current_section == 'elements' and len(parts) >= 4:
                elem_id = int(parts[0])
                node_ids = [int(p) for p in parts[2:] if int(p) != 0]
                elements[elem_id] = [n - 1 for n in node_ids]
            elif current_section == 'displacements' and len(parts) >= 3:
                node_id = int(parts[0])
                displacements[node_id] = [float(parts[1]), float(parts[2]), 0.0]
            elif current_section == 'campo' and current_bloco_nome and len(parts) >= 4:
                elem_id = int(parts[0])
                valores = []
                # Handle potential '*' in values, converting them to 0.0
                for val_str in parts[1:]:
                    if '*' in val_str:
                        valores.append(0.0)
                    else:
                        try:
                            valores.append(float(val_str))
                        except ValueError:
                            # If conversion fails, skip or handle as needed
                            continue
                # Ensure we only take the first 3 values for XX, YY, XY
                if len(valores) >= 3:
                    dados_por_bloco[current_bloco_nome][elem_id] = valores[:3]

    return nodes, elements, displacements, dados_por_bloco


def construir_grid_pyvista(nodes, elements):
    sorted_node_ids = sorted(nodes)
    points = np.array([nodes[i] for i in sorted_node_ids])
    cells = []
    cell_types = []
    for eid in sorted(elements):
        nids = elements[eid]
        if len(nids) >= 3:
            cells.append(len(nids))
            cells.extend(nids)
            cell_types.append(pv.CellType.POLYGON)
    cells_array = np.array(cells, dtype=np.int64)
    cell_types_array = np.array(cell_types, dtype=np.uint8)
    return pv.UnstructuredGrid(cells_array, cell_types_array, points)


def plotar_resultado(filepath, tipo_campo, bloco, componente, escala_deformacao=1.0):
    nodes, elements, displacements, dados_por_bloco = ler_blocos_resultados(filepath)
    grid = construir_grid_pyvista(nodes, elements)

    # Prepara deslocamentos
    disp_array = np.zeros((grid.n_points, 3))
    for nid, val in displacements.items():
        idx = nid - 1
        if 0 <= idx < len(disp_array):
            disp_array[idx] = val

    grid.point_data['Deslocamento_X'] = disp_array[:, 0]
    grid.point_data['Deslocamento_Y'] = disp_array[:, 1]
    grid.point_data['Deslocamento_Magnitude'] = np.linalg.norm(disp_array, axis=1)

    # Aplica deformação
    pontos_deformados = grid.points + disp_array * escala_deformacao
    grid_deformado = grid.copy(deep=True)
    grid_deformado.points = pontos_deformados

    campo_plotar = None
    nome_barra = ""
    if tipo_campo.lower() == 'deslocamentos':
        nome = f"Deslocamento_{componente.upper()}"
        if nome in grid_deformado.point_data:
            campo_plotar = grid_deformado.point_data[nome]
            nome_barra = nome + " (m)"
    else:
        chave_bloco = f"{tipo_campo}_{bloco}"
        if chave_bloco in dados_por_bloco:
            comp_idx = {'XX': 0, 'YY': 1, 'XY': 2}.get(componente.upper(), None)
            if comp_idx is None:
                print(f"Componente inválida: {componente}")
                return
            dados = dados_por_bloco[chave_bloco]
            arr = np.zeros(grid_deformado.n_cells)
            # Mapeia os IDs dos elementos para os índices do array
            elem_id_to_idx = {elem_id: i for i, elem_id in enumerate(sorted(elements))}

            for eid, vals in dados.items():
                if eid in elem_id_to_idx:
                    arr[elem_id_to_idx[eid]] = vals[comp_idx]
            
            grid_deformado.cell_data[f"{tipo_campo}_{componente}_{bloco}"] = arr
            campo_plotar = arr
            nome_barra = f"{tipo_campo}_{componente} ({bloco})"
        else:
            print(f"Bloco não encontrado: {chave_bloco}")
            return

    # Plotagem
    plotter = pv.Plotter()
    plotter.add_mesh(grid_deformado, scalars=campo_plotar, cmap='jet',
                     scalar_bar_args={'title': nome_barra})
    plotter.add_mesh(grid, style='wireframe', color='black', line_width=1)
    plotter.add_axes()
    plotter.add_title(f"Campo: {nome_barra} | Escala deformação: {escala_deformacao}x")
    plotter.show()


# =======================
# CONFIGURAÇÕES DO USUÁRIO
# =======================

filepath = "resultados.txt"
tipo_campo = "Deslocamentos"      # "Tensoes", "Deformacoes", ou "Deslocamentos"
bloco = "Centros dos Elementos"                  # "Centros dos Elementos", "GP 1", ..., "GP 9"
componente = "X"              # "XX", "YY", "XY", ou "X", "Y", "Magnitude" (para deslocamento)
escala_deformacao = 1.0        # Escala visual da deformação

# === EXECUTA ===
plotar_resultado(filepath, tipo_campo, bloco, componente, escala_deformacao)


