import re
from collections import defaultdict

'Altera o formado das tensões e deformaçõe nos pontos de gauss para o modelo ao'
'qual o arquivo plot conegue ler'

def unificar_blocos_gp(linhas):
    blocos = []
    temp = ''
    for linha in linhas:
        if 'GP' in linha and 'Sigma=' in linha:
            temp = linha.strip()
        elif temp:
            temp += ' ' + linha.strip()
            if ')' in linha and 'Epsilon=' in temp:
                blocos.append(temp)
                temp = ''
    return blocos

def extrair_por_gp(blocos_gp):
    dados = defaultdict(list)
    elemento_atual = None
    elem_ids = []
    
    for linha in linhas:
        if '# Elemento Hexa ID:' in linha:
            match = re.search(r'ID:\s+(\d+)', linha)
            if match:
                elemento_atual = int(match.group(1))
                elem_ids.append(elemento_atual)
    
    elemento_iter = iter(elem_ids)
    elemento_id = next(elemento_iter, None)
    contador = 0

    for bloco in blocos_gp:

        match = re.search(
            r'GP\s*(\d+):\s*Sigma=\(([-\d.eE+ ]+)\)\s*Epsilon=\(([-\d.eE+ ]+)\)', bloco
        )
        if match:
            gp = int(match.group(1))
            sigma_raw = match.group(2)
            epsilon_raw = match.group(3)
            try:
                sigma = list(map(float, re.findall(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?', sigma_raw)))
                epsilon = list(map(float, re.findall(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?', epsilon_raw)))
                if len(sigma) == 3 and len(epsilon) == 3:
                    dados[gp].append({
                        'Elemento': elemento_id,
                        'Sigma': sigma,
                        'Epsilon': epsilon
                    })
                    contador += 1
                    if contador % 9 == 0:
                        elemento_id = next(elemento_iter, elemento_id)
            except ValueError as e:
                print(f"Erro de conversão: {e}")
    return dados

with open('resultados.txt', 'r', encoding='utf-8') as f:
    linhas = f.readlines()

blocos_gp = unificar_blocos_gp(linhas)
dados_por_gp = extrair_por_gp(blocos_gp)

with open('por_gp_formatado.txt', 'w', encoding='utf-8') as out:
    for gp in sorted(dados_por_gp.keys()):
        out.write(f"# Tensoes no GP {gp} (Elemento, Sigma_XX, Sigma_YY, Sigma_XY)\n")
        for d in dados_por_gp[gp]:
            s = d['Sigma']
            out.write(f"{d['Elemento']:>5} {s[0]:>25.20f} {s[1]:>25.20f} {s[2]:>25.20f}\n")
        out.write(f"\n# Deformacoes no GP {gp} (Elemento, Epsilon_XX, Epsilon_YY, Epsilon_XY)\n")
        for d in dados_por_gp[gp]:
            e = d['Epsilon']
            out.write(f"{d['Elemento']:>5} {e[0]:>25.20f} {e[1]:>25.20f} {e[2]:>25.20f}\n")
        out.write('\n')

print("✅ Agora foi! O arquivo foi salvo com os dados separados por ponto de Gauss.")
