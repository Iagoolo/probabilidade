
from lightkurve import search_targetpixelfile
import matplotlib.pyplot as plt
import numpy as np

pixelfile = search_targetpixelfile('KIC 8462852', quarter=16).download(quality_bitmask='hardest')
pixelfile.plot(frame=1000)
plt.savefig('video.png')

lc = pixelfile.to_lightcurve(aperture_mask='all')
lc.plot()
plt.savefig('grafico_inicial.png') #tem uma curva nao significa que tem planeta

pixelFile = search_targetpixelfile('KIC 6922244', quarter=4).download()
lc_two = pixelFile.to_lightcurve(aperture_mask=pixelFile.pipeline_mask)
lc_two.plot()
plt.savefig('segunda_grafico.png')

flat_lc = lc_two.flatten(window_length=401)
flat_lc.plot()
plt.savefig('grafico_com_flatten.png')

periodogram = flat_lc.to_periodogram(method='bls', period=np.arange(3, 5, 0.001))
periodogram.plot()
plt.savefig('periodogram.png')

folder_lc = flat_lc.fold(period=periodogram.period_at_max_power)
folder_lc.plot()
plt.savefig('folder.png')

search_result = search_targetpixelfile('Pi Mensae', mission='TESS', sector=1)
tpf = search_result.download(quality_bitmask='default')
lc_pi = tpf.to_lightcurve()
# lc_pi.plot()
# plt.savefig('PiMensae.png')

# lc_pi.scatter()
# plt.savefig('scatter.png')

# aperture_mask = tpf.create_threshold_mask(threshold=10)
# lc = tpf.to_lightcurve(aperture_mask=aperture_mask)
# lc.scatter()
# plt.savefig('scatter_maybe_better.png')

flat_lc = lc_pi.flatten(window_length=1001)
# flat_lc.errorbar()
# plt.savefig('Pierrorbar.png')

mask = (flat_lc.time.value < 1346) | (flat_lc.time.value > 1350)
masked_lc = flat_lc[mask]
masked_lc.errorbar()
plt.savefig('Pierrorbar_better.png')

clipped_lc = masked_lc.remove_outliers(sigma=6)
clipped_lc.errorbar()

# usar o periodograma com bls com o clipped_lc para ver
# aplicar o melhor periodo(fold)
# usar o bin e o binsize para melhor visualizacao
# t0 para mover ao centro


