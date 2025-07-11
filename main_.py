
from lightkurve import search_lightcurve
import matplotlib.pyplot as plt
import numpy as np

# 1. Buscar dados da estrela (exemplo: TOI 700)
lc = search_lightcurve("TOI 700", mission='TESS', sector=1).download()
lc.plot(title="Curva de luz bruta")
plt.savefig("curva_de_luz_b.png")

# 2. Limpar e suavizar a curva
lc_clean = lc.remove_nans().flatten(window_length=401)
lc_clean.plot(title="Curva de luz suavizada")
plt.savefig("curva_luz_s.png")

# 3. Detectar possíveis trânsitos com Box Least Squares
from lightkurve.periodogram import BoxLeastSquares

periods = np.linspace(0.5, 10, 10000)
bls = lc_clean.to_periodogram(method="bls", period=periods)
bls.plot(title="Periodograma BLS")
plt.savefig("periodograma.png")

# 4. Imprimir período mais provável (possível planeta)
print(f"Período mais provável: {bls.period_at_max_power:.4f} dias")


