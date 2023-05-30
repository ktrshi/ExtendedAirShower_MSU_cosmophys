import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from pathlib import Path
import matplotlib as mpl
from iminuit import Minuit
from scipy import integrate
from lets_plot import *
import seaborn as sns
import mplcatppuccin
from mplcatppuccin.palette import load_color

LetsPlot.setup_html()
mpl.style.use("latte")
p0, p1, p2, p3, x0, y0 = 0, 0, 0, 0, 0, 0


def graph(params, transparency=0.1, save=False):
    fig = plt.figure(figsize=(10, 5), dpi=200)
    match params['func']:
        case 'func':
            fig.suptitle(r'$\frac{p_{0}}{1+p_{1}*r+p_{2}*r^{2}}$')
        case 'func_':
            fig.suptitle(r'$\frac{p_{0}}{1+p_{1}*r+p_{2}*r^{1.5}}$')
        case 'func1':
            fig.suptitle(r'$\frac{p_{0}}{1+p_{1}*r+p_{2}*r^{2}+p_{3}*r^{3}}$')

    gs = mpl.gridspec.GridSpec(1, 2)

    ax = fig.add_subplot(gs[0, 0], projection='3d')
    ax.scatter(x_mesh, y_mesh, np.log10(data), color=load_color("mocha", "green"), alpha=transparency)
    ax.plot_surface(x_mesh, y_mesh, np.log10(params['int_func']), cmap='jet')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel(r'$log_{10}$Z')
    ax.set_title('normalised FCN {:.1e}'.format(params['fcn']))

    ax_ = fig.add_subplot(gs[0, 1])
    ax_.scatter(params['rx'], params['diff'], c=np.log10(params['Q']), cmap='jet')
    ax_.set_xlabel('R, m', fontsize=15)
    ax_.set_ylabel('diff / I', fontsize=15)
    ax_.set_title('Relative variation', fontsize=15)

    fig.tight_layout()

    fig.savefig(f'{data_name}{random.randint(1, 1000)}.jpg', dpi=200, bbox_inches='tight')
    # plt.show()


def minim(idx, flag='func'):
    global p0, p1, p2, p3, x0, y0, I1, I2, I3, I4, I5, s
    match flag:
        case 'func':
            m_ = Minuit(fcn=error, p0=10000, p1=10, p2=1, x0=0, y0=0)
            m_.limits = [(0, None), (0, None), (0, None), (-100., 100.), (-100., 100.)]
        case 'func_':
            m_ = Minuit(fcn=error_, p0=10000, p1=10, p2=1, x0=0, y0=0)
            m_.limits = [(0, None), (0, None), (0, None), (-100., 100.), (-100., 100.)]
        case 'func1':
            m_ = Minuit(fcn=error1, p0=5e6, p1=0.2, p2=1, p3=0.1, x0=0, y0=0)
            m_.limits = [(0, None), (0, None), (0, None), (0, None), (-400., 400.), (-400., 400.)]
    # print(m_.migrad())
    m_.migrad()
    print(m_.valid)
    fcn_ = m_.fval
    match flag:
        case 'func' | 'func_':
            p0, p1, p2, x0, y0 = m_.values
        case 'func1':
            p0, p1, p2, p3, x0, y0 = m_.values
    r_ = np.array(list(zip(x_mesh, y_mesh)))
    rx_ = x_mesh.reshape(-1)
    ry_ = y_mesh.reshape(-1)
    match flag:
        case 'func':
            Q_ = func(p0, p1, p2, rx_, ry_, x0, y0)
            int_func_ = func(p0, p1, p2, x_mesh, y_mesh, x0, y0)
        case 'func_':
            Q_ = func_(p0, p1, p2, rx_, ry_, x0, y0)
            int_func_ = func_(p0, p1, p2, x_mesh, y_mesh, x0, y0)
        case 'func1':
            Q_ = func1(p0, p1, p2, p3, rx_, ry_, x0, y0)
            int_func_ = func1(p0, p1, p2, p3, x_mesh, y_mesh, x0, y0)
    diff_ = (Q_ - data.reshape(-1)) / Q_
    sigm2_ = ((Q_ - data.reshape(-1)) ** 2).sum()
    sum = data.sum()
    sum_ = int_func_.sum()
    med = np.median(data[data > 0])
    data[data == 0] = med
    I1[idx] = integrate.quad(f3_, 0, 100, limit=100)[0]
    I2[idx] = integrate.quad(f3_, 0, 200, limit=100)[0]
    I3[idx] = integrate.quad(f3_, 0, 300, limit=100)[0]
    I4[idx] = integrate.quad(f3_, 0, 400, limit=100)[0]
    I5[idx] = integrate.quad(f3_, 0, 500, limit=100)[0]
    ru = pd.DataFrame({'u': data.reshape(-1), 'r': (rx_ ** 2 + ry_ ** 2) ** (1. / 2.)})
    s[idx] = ru[ru['r'] < 100].sum()['u']
    return {'m': m_, 'r': r_, 'rx': rx_, 'ry': ry_, 'Q': Q_, 'int_func': int_func_, 'diff': diff_, 'sigm2': sigm2_,
            'fcn': fcn_, 'func': flag}


def func(p0, p1, p2, x, y, x0, y0):
    r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
    return p0 / (1 + p1 * r + p2 * r ** 2)


def error(p0, p1, p2, x0, y0):
    return ((func(p0=p0, p1=p1, p2=p2, x=x_mesh, y=y_mesh, x0=x0, y0=y0) - data) ** 2).sum()


def func_(p0, p1, p2, x, y, x0, y0):
    r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
    return p0 / (1 + p1 * r + p2 * r ** 1.5)


def error_(p0, p1, p2, x0, y0):
    return ((func_(p0=p0, p1=p1, p2=p2, x=x_mesh, y=y_mesh, x0=x0, y0=y0) - data) ** 2).sum()


def func1(p0, p1, p2, p3, x, y, x0, y0):
    r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
    return p0 / (1 + p1 * r + p2 * r ** 2 + p3 * r ** 3)


def error1(p0, p1, p2, p3, x0, y0):
    return np.sqrt((func1(p0=p0, p1=p1, p2=p2, p3=p3, x=x_mesh, y=y_mesh, x0=x0, y0=y0) - data) ** 2).sum() / np.abs(
        data).sum()


def f3(x, y):
    r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
    return p0 / (1 + p1 * r + p2 * r ** 2 + p3 * r ** 3)


def f3_(r):
    r0 = np.sqrt(x0 ** 2 + y0 ** 2)
    r_t = np.sqrt(abs(r ** 2 - r0 ** 2))
    # r_t = r
    return p0 / (1 + p1 * r_t + p2 * r_t ** 2 + p3 * r_t ** 3)


data_name = '10PeV_p'
I1 = np.zeros(10)
I2 = np.zeros(10)
I3 = np.zeros(10)
I4 = np.zeros(10)
I5 = np.zeros(10)
res = pd.DataFrame(index=np.arange(100), columns=['Id1', 'Id2', 'Id3', 'Id4', 'Id5', 's', 'I1', 'I3', 'del_r'])
s = np.zeros(10)
params = []
table_list = list((Path.cwd() / 'cher_local' / data_name).glob('**/totcnt_*'))
for k, table in enumerate(table_list):
    df = pd.read_csv(table, header=None, sep='\s+')
    df_ref = pd.DataFrame(df.values.reshape(-1, 400).tolist())
    data = df_ref.values
    x_mesh, y_mesh = np.meshgrid(range(-data.shape[0], data.shape[0], 2), range(-data.shape[1], data.shape[1], 2))
    params.append(minim(k, 'func1'))
    # graph(params[k], save=True)
    # plt.show()

# param = minim(9, 'func')
# graph(param, save=True)
# param1 = minim(9, 'func_')
# graph(param1, save=True)
param2 = minim(9, 'func1')
graph(param2, save=True)

df_small = pd.read_csv(Path.cwd() / 'cher_local' / data_name / 'cnfgs', header=0, sep='\s+',
                       names=['1', '2', '3', '4', '5', '6'],
                       encoding='utf8')
events_data = []
for i in range(100):
    data = df_small[['1', '2', '3', '4', '5']][1 + i * 6:6 + i * 6]
    events_data.append([df_small.iloc[i * 6].to_numpy(), data])
for i in range(100):
    data = events_data[i][1].values
    x0, y0 = events_data[i][0][4], events_data[i][0][5]
    p = int(events_data[i][0][1])
    X = np.linspace(-50, 50, 5)
    X_temp = np.linspace(-50, 50, 100)
    Y = np.linspace(-50, 50, 5)
    Y_temp = np.linspace(-50, 50, 100)
    x_mesh, y_mesh = np.meshgrid(X, Y)
    x_mesh_temp, y_mesh_temp = np.meshgrid(X_temp, Y_temp)
    rx = x_mesh.reshape(-1)
    rx_temp = x_mesh_temp.reshape(-1)
    ry = y_mesh.reshape(-1)
    ry_temp = y_mesh_temp.reshape(-1)
    p0, p1, p2, p3, _, _ = params[p - 1]['m'].values
    m = Minuit(error1, p0=p0, p1=p1, p2=p2, p3=p3, x0=x0, y0=y0)
    m.limits = [(0, None), (0, None), (0, None), (0, None), (-50., 50.), (-50., 50.)]
    # print(m.migrad())
    m.migrad()
    p0, p1, p2, p3, x0, y0 = m.values
    Q = func1(p0, p1, p2, p3, rx, ry, x0, y0)
    int_func = func1(p0, p1, p2, p3, x_mesh, y_mesh, x0, y0)
    diff = (Q - data.reshape(-1)) / Q
    diff_temp = None
    fcn = m.fval / m.nfcn
    # m = Minuit(error1, p0=6000, p1=10, p2=1, p3=0.1, x0=events_data[i][0][4], y0=events_data[i][0][5])
    # m.migrad()
    # p0, p1, p2, p3, x0, y0 = m.values
    Id1 = integrate.quad(f3_, 0, 100, limit=100)[0]
    Id2 = integrate.quad(f3_, 0, 200, limit=100)[0]
    Id3 = integrate.quad(f3_, 0, 300, limit=100)[0]
    Id4 = integrate.quad(f3_, 0, 400, limit=100)[0]
    Id5 = integrate.quad(f3_, 0, 500, limit=100)[0]
    # Id6 = integrate.quad(f3, 0, 600)[0]
    # d = abs(I3[p - 1] - Id3) / I3[p - 1]
    # ds = abs(s[p - 1] - Id1) / s[p - 1]
    res['Id1'].iloc[i] = Id1
    res['Id2'].iloc[i] = Id2
    res['Id3'].iloc[i] = Id3
    res['Id4'].iloc[i] = Id4
    res['Id5'].iloc[i] = Id5
    res['s'].iloc[i] = s[p - 1]
    res['I3'].iloc[i] = I3[p - 1]
    res['I1'].iloc[i] = I1[p - 1]
    res['del_r'].iloc[i] = abs(
        ((events_data[i][0][4]) ** 2 + (events_data[i][0][5]) ** 2) ** (1. / 2.) - ((x0) ** 2 + (y0) ** 2) ** (1. / 2.))
    if i % 20 == 0:
        parametrs = {'m': m, 'rx': rx, 'Q': Q, 'int_func': int_func, 'diff': diff, 'fcn': fcn,
                     'func': 'func1'}
        graph(parametrs, transparency=0.5, save=True)
maxim = res['Id3'].max()
mins = res['Id3'].min()
res = res.drop(res.loc[abs(res['Id3']) > maxim].index)
mean = res['Id3'].mean()
rms = math.sqrt((res['Id3'] ** 2).sum() / len(res))
std = res['Id3'].std(ddof=0)
plt.figure()
plt.title(r'$\int 300$, ' + data_name)
bins = np.arange(4e8, 11e8, 2e7)
plt.hist(abs(res['Id3']), bins=bins, alpha=0.5,
         label='Mean = {:.3e}\nRMS = {:.3e}\nStd = {:.3e}'.format(mean, rms, std))
plt.legend()
plt.savefig(f'{data_name}_N300.jpg', dpi=200, bbox_inches='tight')
mean = res['del_r'].mean()
rms = math.sqrt((res['del_r'] ** 2).sum() / len(res))
std = res['del_r'].std(ddof=0)
plt.figure()
plt.title(f'Core location error distribution, {data_name}')
bins = np.arange(0, 7, 0.2)
plt.hist(abs(res['del_r']), bins=bins, alpha=0.5,
         label='Mean = {:.3}\nRMS = {:.3}\nStd = {:.3}'.format(mean, rms, std))
plt.legend()
plt.savefig(f'{data_name}_dR.jpg', dpi=200, bbox_inches='tight')
# plt.show()
print(0)
