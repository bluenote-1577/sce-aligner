import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns
import numpy as np
from scipy import stats
import sys
import pysam

def func(x,a,b,c):
    return 1 - a  / (np.power(x,b))

cmap = sns.color_palette("muted")
cmap = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.style.use('science')

passing_file = "./passing_reads.txt"

passing_reads = set()
f = open(passing_file, 'r')
for line in f:
    spl = line.split()
    div = float(spl[1])
    read = spl[0]
    if div > 0.045 and div < 0.055:
        passing_reads.add(read)

end = 10000


labels = ["SARS‑CoV‑2\n(0.03 Mb)", "E. coli\n(4.64 Mb)", "M. oryzae\n(40.98 Mb)", "D. melanogaster\n(143.73 Mb)", "H. sapiens\n(3054.82 Mb)"]
time_files = ["./cov_results.txt", "./ecoli_results.txt", "./rice_fungus_results.txt", "./fly_results.txt", "./human_results.txt"]
#time_files = ["./ecoli_results.txt"]
#labels = ["SARS‑CoV‑2 (0.03 Mb)", "E. coli (4.64 Mb)", "M. oryzae (40.98 Mb)", "D. melanogaster (143.73 Mb)"]
#time_files = ["./cov_results.txt", "./ecoli_results.txt", "./rice_fungus_results.txt", "./fly_results.txt"]



plt.rcParams.update({'font.size': 9})
fig = plt.figure(figsize=(7, 2))
for i in range(len(time_files)):
    large_gap_c = 0
    anchor_mult_c = 0
    chain_fail_c = 0
    lap_fail = 0
    robust_x = []
    robust_y = [] 

    time_file = time_files[i]
    f = open(time_file, 'r')
    covered = []
    read_lengths = []
    passing_x = []
    passing_y = []
    for line in f:
        spl = line.split()
        name = spl[4].rstrip()
        if name in passing_reads:
            #outlier
            read_lengths.append(float(spl[2]))
            covered.append(float(spl[3]))
            t = float(spl[0])
            if t == -1:
                lap_fail += 1
            if t < -1.5:
                if t == -2:
                    large_gap_c += 1
                elif t == -3:
                    anchor_mult_c += 1

            elif (float(spl[2]) < 1000 or float(spl[3]) > 0.25) and float(spl[3]) > 0.02:
                robust_x.append(float(spl[2]))
                robust_y.append(float(spl[3]))
                passing_x.append(float(spl[2]))
                passing_y.append(float(spl[3]))

            else:
                passing_x.append(float(spl[2]))
                passing_y.append(float(spl[3]))

                chain_fail_c += 1
#                print(name, spl[3])


    x = read_lengths

    print("large_gap_c", large_gap_c, "anchor_mult_c", anchor_mult_c, "chain_fail_c", chain_fail_c, 'lap_fail', lap_fail)
    popt, pcov = curve_fit(func, robust_x, robust_y)
    y = passing_y;
    y_fit = func(passing_x, popt[0], popt[1], popt[2])
    # residual sum of squares
    ss_res = np.sum((y - y_fit) ** 2)

    # total sum of squares
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    # r-squared
    r2 = 1 - (ss_res / ss_tot)
    print(r2)

    n = np.linspace(0,end,1000)
    print(popt, len(robust_x), len(covered))
    ax = plt.subplot(1, len(time_files), i+1)
    plt.plot(passing_x,passing_y, 'o', c = cmap[i], alpha=0.15)
    a = popt[0]
    b = popt[1]
    #string = r'1 - ${a_v}/{x^{b_v}}$'.replace('a_v',str(round(a,3))).replace('b_v',str(round(b,3)))
    string = f'a = {a.round(1)}\nb = {b.round(3)}'
    r2_str = r'$R^2 = val$'.replace('val', str(r2.round(3)))
    plt.plot(n, func(n,popt[0], popt[1],popt[2]), label = string)
    plt.annotate(r2_str, (end*1.5/4.9, 0.4), fontsize = 8)
    plt.xlim([0,end])
    plt.ylim([0,1])
    if i == 2:
        plt.xlabel("Read length (bp)")
    if i == 0:
        plt.ylabel("Chain length / Read length")
    if i != 0:
        plt.yticks([])
    plt.title(labels[i])
    lgnd = plt.legend(fontsize = 8, loc = 'lower center')
#plt.tight_layout()
plt.savefig('./figures/nanopore_recov.pdf',bbox_inches='tight')
plt.show()



