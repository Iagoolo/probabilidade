# Trabalha iago

## Tipos de ruido a serem tratados

### Máscara de abertura (aperture mask)

Na fotometria de séries temporais de estrelas, cada “pixel file” contém um conjunto de pixels ao redor do alvo. 
Nem todos contribuem com sinal — muitos só trazem ruído de fundo ou de campo. 
A máscara de abertura define quais pixels usar para somar a intensidade da estrela em cada quadro.

Por que usar?

- Evita incluir pixels de fundo (céu, estrelas vizinhas) que só aumentam o ruído.
- Melhora  relação sinal-ruído (S/N).

Como gerar?

Threshold mask: inclui pixels com fluxo médio acima de um limiar (e.g. 3σ acima do fundo).

```python
mask = tpf.create_threshold_mask(threshold=3)
lc = tpf.to_lightcurve(aperture_mask=mask)
```

- Pipeline mask: máscara pré-definida pela Lightkurve com base em análises anteriores.
- Manual mask: você desenha uma forma (círculo, polígono) selecionando pixels “à mão”.

### remove_outliers

É uma função do Lightkurve que elimina pontos isolados (outliers) na curva de luz, baseando-se em sigma‐clipping:

```text
lc_clean = lc.remove_outliers(sigma=5, niters=3)
O que faz?
```

- Calcula a média local (ou um modelo) da curva
- Remove pontos que se desviam mais que sigma vezes o desvio-padrão
- Repete até niters iterações.

Por que usar?
• Elimina “spikes” (rupturas abruptas) causados por artefatos de leitura do detector ou eventos momentâneos.

### Raios cósmicos

São partículas de altíssima energia (principalmente prótons) vindas do espaço. Quando atingem o detector de um telescópio:

Efeito no pixel: geram um pulso de carga elétrica MUITO acima do fluxo estelar normal, criando um “spike” no dado.

Como identificar?

- Aparecem como pontos isolados muito acima da curva geral.
- São eliminados pelo sigma‐clipping ou filtros medianos.

Em resumo, são um tipo de ruído pontual que distorce a curva de luz se não for removido.

### Vetores de co-correlação (Co-trending Basis Vectors – CBVs)

São componentes de tendência sistemática extraídos de um grande conjunto de curvas de luz da mesma missão (Kepler/TESS):

O que são?

- Vetores que descrevem variações instrumentais comuns (drifts térmicos, vibrações, variações de foco).
- Cada CBV é um padrão temporal (“modo”) que se repete em muitas estrelas.

Como usar?

```text
corrected_lc = lc.correct(cbv_components=[1,2])
```

- A função ajusta linearmente esses vetores ao seu dado e subtrai as componentes instrumentais.

Por que é útil?

- Remove tendências de longo prazo sem “achatar” o trânsito real.
- Mais sofisticado que o flatten(), pois aproveita informações de toda a missão.

## 1. Escalonamento do Fluxo

O **fluxo** (flux) é a medida de quanta luz a estrela fornece em cada instante, normalmente em unidades de contagens por segundo (cts/s) ou elétrons por segundo (e–/s). Esses valores dependem das características do instrumento e não têm um “zero” ou escala padronizada.

**Escalonar o fluxo** significa dividir todos os valores de fluxo por uma referência (geralmente a mediana ou a média do próprio fluxo), para que a curva resultante fique em uma escala relativa.

- **Fórmula (mediana):**  
  $
    F_{\text{escalonado}}(t) = \frac{F(t)}{\mathrm{mediana}(F)}
  $
- **Fórmula (média):**  
  $
    F_{\text{escalonado}}(t) = \frac{F(t)}{\langle F \rangle}
  $

**Exemplo em Lightkurve (dividindo pela mediana):**

```python
lc_norm = lc.normalize()
```

**Exemplo manual:**

```python
import numpy as np

flux_med = np.median(lc.flux)
lc_norm = lc.copy()
lc_norm.flux     = lc.flux     / flux_med
lc_norm.flux_err = lc.flux_err / flux_med
```

---

## 2. Curva de Luz Adimensional

Uma curva de luz é chamada de **adimensional** quando seus valores de fluxo não têm unidades físicas explícitas, mas sim uma escala relativa sem dimensão (números puros).

- **Por que “adimensional”?**  
  Após o escalonamento, o fluxo é expresso em “vezes a referência” (por exemplo, “vezes a mediana”), portanto não é mais cts/s nem e–/s, mas números sem unidade.  
- **Vantagens:**  
  - Facilita comparar curvas de diferentes estrelas (todas têm base em 1.0).  
  - Destaca automaticamente as quedas de trânsito como desníveis relativos (por exemplo, 0.998 corresponde a –0.2 %).

---

## 3. Medir Profundidade de Trânsito

A **profundidade de trânsito** (δ) é a fração de luz bloqueada pelo planeta quando ele passa na frente da estrela. Em curva adimensional escalonada pela mediana:

$
  δ = 1 - F_{\min}
$

onde $(F_{\min})$ é o valor mínimo do fluxo durante o trânsito.

- **Exemplo numérico:**  
  Se o nível de base é 1.0 e, durante o trânsito, o fluxo cai para 0.995, então  
  $
    δ = 1 - 0.995 = 0.005 = 0.5\%
  $

- **Como extrair em Python:**

```python
import numpy as np

# Supondo lc_folded já dobrada e normalizada:
flux_tr = lc_folded.flux
depth   = 1 - np.min(flux_tr)
print(f"Profundidade de trânsito: {depth*100:.2f}%")
```

- **Relação com o raio do planeta:**  
  $
    δ ≈ \left(\frac{R_p}{R_\star}\right)^2
    \quad\Longrightarrow\quad
    \frac{R_p}{R_\star} = \sqrt{δ}
  $

---

## 4. Recapitulando

1. **Escalonamento do fluxo**  
   Divide o fluxo original pela mediana (ou média), gerando valores relativos.

2. **Curva adimensional**  
   Resultado do escalonamento — números sem unidade, baseados em “1.0”.

3. **Medir profundidade**  
   Diferença entre o nível de base (1.0) e o mínimo do trânsito; representa a fração de luz bloqueada.

