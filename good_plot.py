import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns
import numpy as np
from scipy import stats
import sys
import pysam

def ext_f(x,a,b):
    return a * np.log(x) * np.power(x, 0.08) + b

combined_ext_times = True

cmap = sns.color_palette("muted")
cmap = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.style.use('science')

fs= 16

labels = ["SARS‑CoV‑2\n(0.03 Mb)", "E. coli\n(4.64 Mb)", "M. oryzae\n(40.98 Mb)", "D. melanogaster\n(143.73 Mb)", "H.sapiens\n(3054.82 Mb)"]
labels_nobp = ["SARS‑CoV‑2 (k = 16)", "E. coli (k = 24)", "M. oryzae (k = 27)", "D. melanogaster (k = 29)", "H. sapiens (k = 34)"]
time_files = ["./cov_results.txt", "./ecoli_results.txt", "./rice_fungus_results.txt", "./fly_results.txt", "./human_results.txt"]
lengths = [30000, 4.64 * 10e6, 40.98 * 10e6, 143.73 * 10e6 , 3054.82 * 10e6]
passing_file = "./passing_reads.txt"

ends = [24000, 90e12]

passing_reads = set()
f = open(passing_file, 'r')
for line in f:
    spl = line.split()
    div = float(spl[1])
    read = spl[0]
    if div > 0.045 and div < 0.055:
        passing_reads.add(read)

if combined_ext_times:
    end = 75000
    yend = 0.06
    slopes_ext = []
    slopes_chain = []
    plt.figure(figsize=(7, 3.5))
    plt.rcParams.update({'font.size': 9})
    for i in range(len(time_files)):
        time_file = time_files[i]
        f = open(time_file, 'r')
        length = []
        time_ext = []
        time_chain = []
        for line in f:
            spl = line.split()
            name = spl[4].rstrip()
            time = float(spl[0])
            if time < 0:
                continue
            if name in passing_reads:
                length.append(float(spl[2]))
                time_ext.append(float(spl[0]) * 10e3)
                time_chain.append(float(spl[1]) * 10e3)


        xrange = np.array(list(range(0,100000)));
        #print(stats.linregress(length,time), 'pearson')
        (slope,intercept) = stats.siegelslopes(time_ext,length)
        print(slope, time_file)
        slopes_ext.append(slope)
        ax = plt.subplot(2, len(time_files), i+1)
        print(len(time_ext))
        plt.plot(length, time_ext, 'o', c = cmap[i], markersize = 3, alpha= 0.15)
        plt.plot(xrange, intercept+slopes_ext[i] * np.array(xrange), 'r-', c = cmap[i])
        plt.title(labels[i])
        if i ==0:
            plt.ylabel("Extension time (ms)")
            plt.xlim([0,end/60])
            plt.ylim([0,yend/60 * 10e3])
        else:
            plt.xlim([0,end])
        #plt.xlim([0,np.min([np.max(length),25000])])
            plt.ylim([0,yend* 10e3])
#        plt.xlim([0,15000])

        (slope,intercept) = stats.siegelslopes(time_chain,length)
        print(slope, time_file)
        slopes_chain.append(slope)
        ax = plt.subplot(2, len(time_files), i+1+ len(time_files))
        plt.plot(length, time_chain, 'o', c = cmap[i], markersize = 3, alpha = 0.1)
        plt.plot(xrange, intercept + slope * np.array(xrange), 'r-', c = cmap[i])
        if i == 2:
            plt.xlabel("Read length (bp)")
        if i == 0:
            plt.ylabel("Chaining time (ms)")
            plt.xlim([0,end/60])
            plt.ylim([0,yend/60/10 * 10e3])
        else:
        #plt.xlim([0,np.min([np.max(length),25000])])
            plt.xlim([0,end])
            plt.ylim([0,yend/10* 10e3])
#        plt.xlim([0,15000])

    plt.tight_layout()
    plt.savefig("figures/combined_ext_plot.pdf")
    plt.show()

    plt.rcParams.update({'font.size': 9})
    fig = plt.figure(figsize=(4, 3))


    popt, pcov = curve_fit(ext_f, lengths, slopes_ext)

    y = slopes_ext
    y_fit = ext_f(lengths, popt[0], popt[1])
    # residual sum of squares
    ss_res = np.sum((y - y_fit) ** 2)
    # total sum of squares
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    # r-squared
    r2 = 1 - (ss_res / ss_tot)
    print(r2)

    res = stats.linregress(np.log(lengths), slopes_ext)
    print(res)

    sp = np.logspace(np.log(ends[0]),np.log(ends[1]),num=1000, base = 2.713,endpoint=True)

    #plt.plot(ends, res.slope * np.log(np.array(ends)) + res.intercept, '--', c = 'black', label = f"y = {res.slope:.2E} log(x) - {np.abs(res.intercept):.2E}")
    #plt.plot(ends, res.slope * np.log(np.array(ends)) + res.intercept, '--', c = 'black', label = f"y = c log(x) x^{C \alpha} + d")

    l1 = r"$y = A_1 \log(x) + B_1$" + '\n' +  r"$R^2 = rval$".replace('rval', str(round(res.rvalue**2,3)))
    l2 = r"$y = A_2 \log(x) x^{C \alpha}  + B_2$" + '\n'+ "$R^2 = rval$".replace('rval',str(round(r2,3)))
    plt.plot(ends, res.slope * np.log(np.array(ends)) + res.intercept, linestyle='dashed', c = 'black', label = l1)
    plt.plot(sp, ext_f(sp,popt[0],popt[1]), linestyle='dotted', c = 'black', label = l2)

    for i in range(len(time_files)):
        plt.plot(lengths[i],slopes_ext[i], 'o', c = cmap[i])
        if i == 1:
            plt.annotate(labels_nobp[i], (lengths[i]*1.4, slopes_ext[i]*0.9))
        else:
            plt.annotate(labels_nobp[i], (lengths[i]*1.4, slopes_ext[i]))
    plt.legend()
    plt.xscale('log')
    plt.xlim(ends)
    plt.xlabel("Genome length (bp)")
    plt.ylabel("Slope of regression (extension runtime)")
    plt.ylim([0,0.012])
    plt.savefig("figures/log_coeff_plot_ext.pdf")
    plt.show()

    #fig = plt.figure(figsize=(8, 6))
    #ends = [28000, 10e11]
    #res = stats.linregress(np.log(lengths), slopes_chain)
    #plt.plot(ends, res.slope * np.log(np.array(ends)) + res.intercept, '--', c = 'black', label = f"y = {res.slope:.2E} log(x) + {res.intercept:.2E}")
    #for i in range(len(time_files)):
    #    plt.plot(lengths[i],slopes_chain[i], 'o', c = cmap[i])
    #    plt.annotate(labels_nobp[i], (lengths[i]*1.4, slopes_chain[i]))
    #plt.legend()
    #plt.xscale('log')
    #plt.xlim(ends)
    #plt.xlabel("Genome length (bp)")
    #plt.ylabel("Slope of regression (chaining runtime)")
    #plt.savefig("figures/log_coeff_plot_chain.pdf")
    #plt.show()
