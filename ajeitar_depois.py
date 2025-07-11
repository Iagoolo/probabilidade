# Cell 1: Imports e Configurações
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lightkurve import TessTargetPixelFile
from scipy.signal import savgol_filter
from scipy import stats
import batman
import emcee

# Esta célula importa bibliotecas e configura parâmetros gerais.

# Cell 2: Definição de process_lightcurve

def process_lightcurve(
    fname,
    smooth_window=51,
    smooth_poly=2,
    bls_periods=np.arange(3, 5, 0.001)
):
    """
    1) Abre o arquivo .fits do TESS (TessTargetPixelFile)
    2) Extrai light curve usando pipeline_mask
    3) Remove outliers (sigma=5)
    4) Faz flatten padrão
    5) Suaviza sinal
    6) Calcula periodograma BLS

    Retorna dict com:
      - file, time, flux, flux_err, t0, best_period,
        mean_flux, std_flux, n_points
    """
    # 1) Leitura do PixelFile a partir do nome no disco
    pixel_file = TessTargetPixelFile(fname)

    # 2) Extração da lightcurve com máscara do pipeline
    lc = pixel_file.to_lightcurve(aperture_mask=pixel_file.pipeline_mask)

    # 3) Remoção de outliers (fluxos muito fora da média)
    lc = lc.remove_outliers(sigma=5)

    # 4) Flatten padrão (remove tendência)
    flat = lc.flatten()

    # 5) Suavização opcional do fluxo
    flux_smooth = savgol_filter(flat.flux, window_length=smooth_window, polyorder=smooth_poly)

    # 6) Periodograma BLS para detectar trânsitos
    periodogram = flat.to_periodogram(method='bls', period=bls_periods)
    best_p = periodogram.period_at_max_power.value

    # Estatísticas do fluxo suavizado
    mean_flux = np.mean(flux_smooth)
    std_flux  = np.std(flux_smooth, ddof=1)
    n_points  = len(flux_smooth)

    return {
        'file': fname.split('/')[-1],
        'time': flat.time.value,
        'flux': flux_smooth,
        'flux_err': flat.flux_err,
        't0': periodogram.transit_time_at_max_power.value,
        'best_period': best_p,
        'mean_flux': mean_flux,
        'std_flux': std_flux,
        'n_points': n_points
    }

# Cell 3: Definição de fit_transit

def fit_transit(
    time,
    flux,
    flux_err,
    period,
    t0,
    r_init,
    a_init,
    inc_init=89.0,
    n_walkers=32,
    n_steps=5000,
    burn=1000
):
    """
    Ajusta modelo de trânsito com batman e emcee.
    Retorna amostras do posterior para [t0, rp, a, inc].
    """
    params = batman.TransitParams()
    params.t0        = t0
    params.per       = period
    params.rp        = r_init
    params.a         = a_init
    params.inc       = inc_init
    params.ecc       = 0.0
    params.w         = 90.0
    params.limb_dark = "quadratic"
    params.u         = [0.3, 0.2]

    def log_prior(theta):
        t0_, rp_, a_, inc_ = theta
        if 0 < rp_ < 0.5 and 1 < a_ < 100 and 80 < inc_ < 90:
            return 0.0
        return -np.inf

    def log_likelihood(theta, t, f, ferr):
        t0_, rp_, a_, inc_ = theta
        params.t0, params.rp, params.a, params.inc = t0_, rp_, a_, inc_
        model = batman.TransitModel(params, t).light_curve(params)
        return -0.5 * np.sum(((f - model) / ferr)**2)

    def log_posterior(theta, t, f, ferr):
        lp = log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_likelihood(theta, t, f, ferr)

    pos = np.array([t0, r_init, a_init, inc_init])
    pos = pos + 1e-4 * np.random.randn(n_walkers, len(pos))

    sampler = emcee.EnsembleSampler(
        n_walkers, len(pos[0]), log_posterior,
        args=(time, flux, flux_err)
    )
    sampler.run_mcmc(pos, n_steps, progress=True)

    return sampler.get_chain(discard=burn, flat=True)

# Cell 4: Execução do Pipeline e Análise Inicial

def run_pipeline(path_pattern='path/para/seus_fits/*.fits'):
    files = glob.glob(path_pattern)
    summary = []

    for fname in files:
        res = process_lightcurve(fname)
        depth = 1 - np.min(res['flux']) / np.median(res['flux'])
        r_init = np.sqrt(depth)
        a_init = (res['best_period']/365.25)**(2/3) * 215

        samples = fit_transit(
            res['time'], res['flux'], res['flux_err'],
            res['best_period'], res['t0'], r_init, a_init
        )

        med = np.median(samples, axis=0)
        summary.append({
            'file':   res['file'],
            'rp/Rs':  med[1],
            'a/Rs':   med[2],
            'inc':    med[3],
            'period': res['best_period']
        })

    df = pd.DataFrame(summary)
    df.to_csv('catalogo_exoplanetas.csv', index=False)

    print(df.describe())
    plt.figure()
    plt.hist(df['rp/Rs'], bins=20)
    plt.xlabel('Rp/Rs')
    plt.ylabel('Contagem')
    plt.title('Distribuição de Rp/Rs na amostra')
    plt.show()

    return df

# Cell 5: Análise de Correlações e Occurrence Rate

def analyze_results(df):
    import seaborn as sns
    vars_to_corr = ['rp/Rs', 'a/Rs', 'period']
    print("Matriz de correlação:")
    print(df[vars_to_corr].corr())
    sns.pairplot(df[vars_to_corr])
    plt.suptitle("Scatter e distribuições marginais", y=1.02)
    plt.show()

    df['transit_prob'] = 1 / df['a/Rs']
    occurrence_rate = len(df) / df['transit_prob'].sum()
    print(f"Occurrence rate aproximada: {occurrence_rate:.3f} planetas por estrela")

# Uso:
# df = run_pipeline()
# analyze_results(df)
