
from lightkurve import search_lightcurve
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

# coord = SkyCoord(ra=99.083374*u.degree, dec=11.866667*u.degree)
region = search_lightcurve("00h42m44.3s +41d16m9s", radius=10*u.arcmin, mission="TESS")
# region.plot(title="Regiao")
# plt.savefig("todos.png")

cont = 0

for star in region[:500]:
    cont += 1

print(cont)
