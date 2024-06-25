#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import statistics
import random
import datetime
import numpy as np
import pandas as pd
# avoid error 'Invalid DISPLAY variable' by specifying matplotlib backend before importing pyplot
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main():

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
        
    sep1 = "#" * 10
    sep2 = "=" * 10

    ############### collect data ################
    print("{sep} read data {sep}".format(sep=sep1))

    # site = (seqid, position)
    
    print("{sep} read bed_file".format(sep=sep2))
    # cats = list os categories
    # items = list of items
    # mir_pos = list of position (should be [1 .. 21])
    # site2cats = dict, key = site, value = list of cats
    # site2item = dict, key = site, value = list of items
    # site2mir_pos = dict, key = site, value = list pos
    cats, site2cats, items, site2items, mir_poss, site2mir_pos, n_bed_sites = read_bed_file(args.in_bed)
    print("{} categories : {}".format(len(cats), cats))

    print("{sep} read pi_file".format(sep=sep2))
    # all_pi = dict, key = site, value = pi
    all_pi, n_in_pi = read_pi_file(args.in_pi, site2items)

    ################ compute values ################
    print("{sep} Compute observed mean pi {sep}".format(sep=sep1))
    # site_vector = list of all sites [site = (seqid, position)]]
    # pi_vector = list of all pi values
    site_vector = all_pi.keys()
    observed_pi_vector = all_pi.values()
    n_pi_eq0 = sum([1 for i in observed_pi_vector if i == 0])
    n_pi_gt0 = sum([1 for i in observed_pi_vector if i > 0])
    n_pi_lt0 = sum([1 for i in observed_pi_vector if i < 0])
    mean_pi = statistics.mean(observed_pi_vector)
    print("{:,} sites, {:,} pi, {:,} pi = 0, {:,} pi > 0, {:,} pi < 0".format(len(site_vector), len(observed_pi_vector), n_pi_eq0, n_pi_gt0, n_pi_lt0))
    # per cat
    print("{sep} Compute observed mean pi per category {sep}".format(sep=sep2))
    observed_pi_list_per_cat = make_pi_list_per_cat(observed_pi_vector, site_vector, site2cats, cats)
    effectif_per_cat = dict()
    effectif_polymorphe_per_cat = dict()
    observed_mean_pi_per_cat = dict()
    for cat in cats:
        effectif_per_cat[cat] = len(observed_pi_list_per_cat[cat])
        effectif_polymorphe_per_cat[cat] = len([ i for i in observed_pi_list_per_cat[cat] if i > 0])
        if effectif_per_cat[cat] < 3:
            observed_mean_pi_per_cat[cat] = "NA"
        else:
            observed_mean_pi_per_cat[cat] = statistics.mean(observed_pi_list_per_cat[cat])


    # per item
    print("{sep} Compute observed mean pi per item {sep}".format(sep=sep2))
    observed_pi_list_per_item = make_pi_list_per_cat(observed_pi_vector, site_vector, site2items, items)
    effectif_per_item = dict()
    effectif_polymorphe_per_item = dict()
    observed_mean_pi_per_item = dict()
    for item in items:
        effectif_per_item[item] = len(observed_pi_list_per_item[item])
        effectif_polymorphe_per_item[item] = len([ i for i in observed_pi_list_per_item[item] if i > 0])
        if effectif_per_item[item] < 3:
            observed_mean_pi_per_item[item] = "NA"
        else:
            observed_mean_pi_per_item[item] = statistics.mean(observed_pi_list_per_item[item])

    # per mir_pos
    print("{sep} Compute observed mean pi per mir_pos {sep}".format(sep=sep2))
    observed_pi_list_per_mir_pos = make_pi_list_per_cat(observed_pi_vector, site_vector, site2mir_pos, mir_poss)
    effectif_per_mir_pos = dict()
    effectif_polymorphe_per_mir_pos = dict()
    observed_mean_pi_per_mir_pos = dict()
    for mir_pos in mir_poss:
        effectif_per_mir_pos[mir_pos] = len(observed_pi_list_per_mir_pos[mir_pos])
        effectif_polymorphe_per_mir_pos[mir_pos] = len([ i for i in observed_pi_list_per_mir_pos[mir_pos] if i > 0])
        if effectif_per_mir_pos[mir_pos] < 3:
            observed_mean_pi_per_mir_pos[mir_pos] = "NA"
        else:
            observed_mean_pi_per_mir_pos[mir_pos] = statistics.mean(observed_pi_list_per_mir_pos[mir_pos])

    ################ create csv results files ################
    # - date, input file and result summary first lines of result csv files
    header_tab = ["# {}".format(datetime.datetime.today().strftime('%Y-%m-%d-%H:%M:%S')),
                  "bed_file={}".format(args.in_bed),
                  "pi_file={}".format(args.in_pi),
                  "items={:,}".format(len(items)),
                  "cats={:,}".format(len(cats)),
                  "sites={:,}".format(n_bed_sites),
                  "pi_values={:,}/{:,}".format(len(site_vector), n_in_pi),
                  "mean pi={:,}".format(mean_pi),
    ]
    header_line = "\n# ".join(header_tab) + "\n"

    # - observed values per cat
    csv_file = "{}/observed_per_cat.csv".format(args.out_dir)
    print("{sep} write {out_file} {sep}".format(sep=sep2, out_file=csv_file))
    with open(csv_file, 'w') as f_out:
        f_out.write(header_line)
        tab = ["cat", "n_sites", "n_sites_polymorphes", "mean_pi"]
        f_out.write("{}\n".format("\t".join(tab)))
        for cat in cats:
            tab = [cat, effectif_per_cat[cat], effectif_polymorphe_per_cat[cat], observed_mean_pi_per_cat[cat]]
            tab = [str(i) for i in tab]
            f_out.write("{}\n".format("\t".join(tab)))

    # - observed values per item
    csv_file = "{}/observed_per_item.csv".format(args.out_dir)
    print("{sep} write {out_file} {sep}".format(sep=sep2, out_file=csv_file))
    with open(csv_file, 'w') as f_out:
        f_out.write(header_line)
        tab = ["item", "size", "n_sites", "n_sites_polymorphes", "mean_pi"]
        f_out.write("{}\n".format("\t".join(tab)))
        for item in items:
            tmp = item.split('_')
            start = int(tmp[-2])
            end = int(tmp[-1])
            size = end - start + 1
            tab = [item, size, effectif_per_item[item], effectif_polymorphe_per_item[item], observed_mean_pi_per_item[item]]
            tab = [str(i) for i in tab]
            f_out.write("{}\n".format("\t".join(tab)))

    # - observed values per mir_pos
    csv_file = "{}/observed_per_mir_pos.csv".format(args.out_dir)
    print("{sep} write {out_file} {sep}".format(sep=sep2, out_file=csv_file))
    with open(csv_file, 'w') as f_out:
        f_out.write(header_line)
        tab = ["pos", "n_sites", "n_sites_polymorphes", "mean_pi"]
        f_out.write("{}\n".format("\t".join(tab)))
        for mir_pos in mir_poss:
            tab = [mir_pos, effectif_per_mir_pos[mir_pos], effectif_polymorphe_per_mir_pos[mir_pos], observed_mean_pi_per_mir_pos[mir_pos]]
            tab = [str(i) for i in tab]
            f_out.write("{}\n".format("\t".join(tab)))


    ################ simul stuff ################
    if args.simul > 0:
        # compute simulated mean pi
        print("{sep} Compute simulated mean pi per category {sep}".format(sep=sep1))
        simul_mean_pi_per_cat = make_simul_mean_pi_per_cat(args.simul, effectif_per_cat, site_vector, observed_pi_vector, site2cats, cats)

        # create output files for simulated values
        for cat in cats:
            out_file = "{}/simul_{}_{}".format(args.out_dir, args.simul, cat)
            print("{sep} write {out_file} {sep}".format(sep=sep2, out_file=out_file))
            with open(out_file, 'w') as f_out:
                for pi in simul_mean_pi_per_cat[cat]:
                    f_out.write("{}\n".format(pi))
                    
        # create csv results files
        #       - simulated values : n_iter, n_sites, n_poly_sites, n_sites_le_mean, one sided and two-sided confidence interval
        csv_file = "{}/simul_{}.csv".format(args.out_dir, args.simul)
        print("{sep} write {out_file} {sep}".format(sep=sep2, out_file=csv_file))
        percentiles = [2.5, 97.5, 5, 95]
        # count = number of simulated mean pi values less or equal to observed mean pi value
        with open(csv_file, 'w') as f_out:
            f_out.write(header_line)
            tab = ["cat", "n_sites", "n_sites_polymorphes", "observed_mean_pi", "n_iter", "count0", "count", "simulated_mean_pi"]
            for percentile in percentiles:
                tab.append("percentile {}".format(percentile))
            f_out.write("{}\n".format("\t".join(tab)))
            for cat in cats:
                simul_mean = statistics.mean(simul_mean_pi_per_cat[cat])
                thr = observed_mean_pi_per_cat[cat]
                count = sum([1 for pi in simul_mean_pi_per_cat[cat] if pi <= thr])
                count0 = sum([1 for pi in simul_mean_pi_per_cat[cat] if pi == 0])
                tab = [cat, effectif_per_cat[cat], effectif_polymorphe_per_cat[cat], observed_mean_pi_per_cat[cat], args.simul, count0, count, simul_mean]
                for percentile in percentiles:
                    tab.append(np.percentile(simul_mean_pi_per_cat[cat], percentile))
                tab = [str(i) for i in tab]
                f_out.write("{}\n".format("\t".join(tab)))

    ################ plot results to pdf ################
    if args.plot == "Y":
        pdf_file = "{}/simul_{}_plot.pdf".format(args.out_dir, args.simul)
        print("{sep} Plot results to {pdf} {sep}".format(sep=sep1, pdf=pdf_file))
        pdf = PdfPages(pdf_file)
        # limit for x axis
        xlim = (args.xmin, args.xmax)

        print("{sep} observed mean pi values as boxplot {sep}".format(sep=sep2))
        plot_observed_values(pdf, cats, effectif_per_cat, observed_pi_list_per_cat)

        print("{sep} simulated distrib as density curve {sep}".format(sep=sep2))
        plot_simulated_distrib(pdf, cats, effectif_per_cat, simul_mean_pi_per_cat, xlim)

        print("{sep} observed mean pi value + simulated distrib as density curve {sep}".format(sep=sep2))
        for cat in cats:
            plot_observed_value_and_simulated_distrib(pdf, cat, effectif_per_cat[cat], observed_mean_pi_per_cat[cat], simul_mean_pi_per_cat[cat], xlim)

        # finish pdf
        d = pdf.infodict()
        d['Title'] = 'Pi Analysis per annotation'
        d['Author'] = u'XXX'
        d['Subject'] = 'blabla'
        d['Keywords'] = 'foo bar'
        d['CreationDate'] = datetime.datetime.today()
        d['ModDate'] = datetime.datetime.today()
        pdf.close()

########################## utils #############################"

def make_pi_list_per_cat(pi_vector, site_vector, site2cats, cats):
    """ 
    create lists of pi per category
    input : 
    pi_vector = list of pi values
    site_vector = list of sites
    site2cats = dict, key = site, value = list of cats
    output :
    pi_list_per_cat = dict, key = cat, value = list of pi values
    examples : 
    cats['stem'] = (0.1, 0.3, 0.8, ...)
    cats['loop'] = (0.3, 0.7 ...)
    """

    pi_list_per_cat = dict()
    for cat in cats:
        pi_list_per_cat[cat] = list()

    ok = 0
    missing = 0
    for pi, site in zip(pi_vector, site_vector):
        seqid, pos = site
        if site not in site2cats:
            missing += 1
            continue
        ok += 1
        for cat in site2cats[site]:
            pi_list_per_cat[cat].append(pi)

    #if missing > 0:
        #print("mk_pi_list_per_cat : {} missing sites".format(missing))
    if ok == 0:
        print("Error : mk_pi_list_per_cat : 0 sites found")
    return pi_list_per_cat

def make_simul_mean_pi_per_cat(n_iter, effectif_per_cat, site_vector, observed_pi_vector, site2cats, cats):
    """ compute simulated mean pi per cat
    for each iteration : 
        shuffle pi_vector
        split pi values to categories
        compute mean pi per category
    """

    # simul_mean_pi_per_cat = dict(), key = cat, value = list of n_iter mean pi
    simul_mean_pi_per_cat = dict()
    for cat in cats:
        simul_mean_pi_per_cat[cat] = list()
    for i in range(n_iter):
        # first copy list
        simul_pi_vector = list(observed_pi_vector)
        # then shuffle in place
        random.shuffle(simul_pi_vector)
        simul_pi_list_per_cat = make_pi_list_per_cat(simul_pi_vector, site_vector, site2cats, cats)
        for cat in cats:
            pi_list = simul_pi_list_per_cat[cat]
            pi_mean = statistics.mean(pi_list)
            simul_mean_pi_per_cat[cat].append(pi_mean)

    return simul_mean_pi_per_cat

########################## read input files #############################"
def norm_scaffold(name):
    """ 
    normalize scaffold name : 
    keep only the scaffold name
    remove the second part (|sizeXXX in gff or :XXX in VCF)
    """
    return name.split('|')[0].split(':')[0]

def read_pi_file(pi_file, site2items):
    """ Collect all pi from a file
    keep only pi values for sites in bed items
    input : 
    pi file = one line per site, 3 cols : seqid position pi
    output : 
    all_pi = fict, key = site (ie pair [seqid, pos]), value = pi
    """

    n_in = 0
    all_pi = dict()
    with open(pi_file) as f_in:
        for l in f_in:
            if l.startswith('CHROM'):
                continue
            n_in += 1
            seqid, pos, pi = l.strip().split()
            seqid = norm_scaffold(seqid)
            site = (seqid, int(pos))
            if site in site2items:
                all_pi[site] = float(pi)
    print("{:,}/{:,} pi values from {}".format(len(all_pi), n_in, pi_file))
    return all_pi, n_in

def read_bed_file(bed_file):
    """ create lists of sites from bed file
    site = (seqid, position)
    overlap is treated through a list of cat for each site
    input : 
       bed_file [line = seqid, start, end, category]
    output : 
       cats = list of categories
       site2cats = dict, key = site, value = list of categories
    example : 
    input bed line 
    seq1  1   100   control
    seq1  5    30   exon
    seq1  6     8   stem
    seq1  9    12   loop
    output : 
    cats = [control, exon, stem, loop]
    site2cats[(seq1, 1)] = [control]
    site2cats[(seq1, 2)] = [control]
    ...
    site2cats[(seq1, 5)] = [control, exon]
    site2cats[(seq1, 6)] = [control, exon, stem]
    site2cats[(seq1, 7)] = [control, exon, stem]
    site2cats[(seq1, 8)] = [control, exon, stem]
    site2cats[(seq1, 9)] = [control, exon, loop]
    ...
    ...

    items = list of items
    site2item = dict, key = site, value = list of items

    mir_pos = list of position (should be [1 .. 21])
    site2mir_pos = dict, key = site, value = list pos

    check for duplicate items (seqid+start+end) if args --skip_bed_dup == Y
    write warning for bed duplicates to a file
    """

    warning_dup_file = "{}/warning_duplicated_items.txt".format(args.out_dir)
        
    cats = set()
    site2cats = dict()

    items = list()
    site2items = dict()

    mir_poss = set()
    site2mir_pos = dict()
    
    n_sites = 0
    with open(bed_file) as f_in, open(warning_dup_file, 'w') as f_out:
        for l in f_in:
            seqid, start, end, cat = l.strip().split()
            seqid = norm_scaffold(seqid)
            cats.add(cat)
            item = "{}_{}_{}_{}".format(cat, seqid, start, end)
            if item in items:
                # always report this duplication
                if args.skip_bed_dup == "Y" :
                    # skip duplications base on skip_bed_dup
                    f_out.write("{} skipped\n".format(item))
                    continue
                else:
                    f_out.write("{} keep\n".format(item))
            items.append(item)
            for pos in range(int(start), int(end) + 1):
                site = (seqid, pos)
                n_sites += 1
                if site not in site2cats:
                    site2cats[site] = list()
                if site not in site2items:
                    site2items[site] = list()
                site2cats[site].append(cat)
                site2items[site].append(item)
            # work on mir_pos if size in [19-23]
            size = int(end) - int(start) + 1
            #print("{} {}".format(item, size))
            if size < 19:
                continue
            if size > 24:
                continue
            mir_pos = 0
            for pos in range(int(start), int(end) + 1):
                mir_pos += 1
                mir_poss.add(mir_pos)
                site = (seqid, pos)
                if site not in site2mir_pos:
                    site2mir_pos[site] = list()
                site2mir_pos[site].append(mir_pos)
    
    cats = sorted(list(cats))
    items = sorted(list(items))
    mir_poss = sorted(list(mir_poss))
    
    print("{:,} items, {:,} sites, {:,} cat, {:,} mir_pos from file {}".format(len(items), n_sites, len(cats), len(mir_poss), bed_file))
    return(cats, site2cats, items, site2items, mir_poss, site2mir_pos, n_sites)


def plot_observed_values(pdf, cats, effectif_per_cat, observed_pi_list_per_cat):
    """ one page for observed values as boxplot """
    
    data = [ observed_pi_list_per_cat[cat] for cat in cats]
    labels = ["{} ({:,})".format(cat, effectif_per_cat[cat]) for cat in cats]
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    bp = ax.boxplot(data)
    ax.set_xticklabels(labels, rotation=45, fontsize=6)
    plt.suptitle("Observed pi values per category")
    pdf.savefig(fig)

def plot_simulated_distrib(pdf, cats, effectif_per_cat, simul_mean_pi_per_cat, xlim):
    """ one page for simulated values as density curve """
    
    plt.figure()
    # create ax using only first category
    cat = cats[0]
    simulated_mean_pi = statistics.mean(simul_mean_pi_per_cat[cat])
    label = "{} (n={:,}, mean={:.5f})".format(cat, effectif_per_cat[cat], simulated_mean_pi)
    df = pd.DataFrame(simul_mean_pi_per_cat[cat], columns=[label, ])
    ax = df.plot.density(xlim=xlim)
    # use ax for other categories
    for cat in cats[1:]:
        simulated_mean_pi = statistics.mean(simul_mean_pi_per_cat[cat])
        label = "{} (n={:,}, mean={:.5f})".format(cat, effectif_per_cat[cat], simulated_mean_pi)
        df = pd.DataFrame(simul_mean_pi_per_cat[cat], columns=[label, ])
        df.plot.density(ax=ax)
    plt.suptitle("Simulated mean pi values, iterations={:,}".format(args.simul))
    pdf.savefig()
    plt.close()
    
def plot_observed_value_and_simulated_distrib(pdf, cat, effectif, observed_mean_pi, simul_mean_pi_list, xlim,):
    """ one page with simulated values as density curve + observed value as vertical bar """
    
    simulated_mean_pi = statistics.mean(simul_mean_pi_list)
    #print("{} : observed mean pi value + simulated distrib as density curve [{:.5f} / {:.5f}]".format(cat, observed_mean_pi, simulated_mean_pi))

    label_observed = "observed mean pi {:.5f}".format(observed_mean_pi)
    label_simulated = "simulated mean pi {:.5f}".format(simulated_mean_pi)

    # mix matplotlib and pandas since I don't know how to plot density with matplotlib
    plt.figure()
    # 1) simulated mean pi values as density curve
    df = pd.DataFrame(simul_mean_pi_list, columns=['simulated means pi', ])
    ax = df.plot.density(xlim=xlim)
    # mean pi value for simulated mean distribution as vertical trait
    plt.axvline(x=simulated_mean_pi, color='blue', linestyle='--', label=label_simulated)
    # observed mean pi value as vertical trait 
    plt.axvline(x=observed_mean_pi, color='red', linestyle='--', label=label_observed)
    plt.suptitle("Observed and simulated mean pi values for '{}'".format(cat))
    plt.title("n={:,} sites, iterations={:,}".format(effectif, len(simul_mean_pi_list)))
    plt.legend()
    pdf.savefig()
    plt.close()

if __name__ == '__main__':
    description = """ 
    pi values analysis from one bed file and one "pi" file

    Input files : 
    - bed file : 4 columns : CHROM START END CATEGORY
    - pi file : 3 colomns : CHROM POS PI_VALUE

    item = one line in bed file
    category = col 4 in bed file (eg stem, loop, mir, mir_star ...)
    items can overlap since each site has a list of categories
    mir_pos : [1-23] for item with size [19-23] nuc
    
    For observed pi values : 
     - compute mean pi value per category, per item and per mir position
    For each simul iteration :
    - shuffle pi_value vector
    - compute mean pi value per category based on site's categories

    Output files : in out_dir folder
    - observed_per_cat.csv : cat	n_sites	n_sites_polymorphes	mean_pi
    - observed_per_item.csv : item	size	n_sites	n_sites_polymorphes	mean_pi
    - observed_per_mir_pos.csv : pos	n_sites	n_sites_polymorphes	mean_pi
    - simul_N.csv : cat	n_sites	n_sites_polymorphes	observed_mean_pi	n_iter	count0	count	simulated_mean_pi	percentile 2.5	percentile 97.5	percentile 5	percentile 95
    - simul_N_category : N lines with the mean pi value for each simul iteration
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--in_bed", help="input bed file", required=True)
    parser.add_argument("--in_pi", help="input pi file", required=True)
    parser.add_argument("--out_dir", help="output dir for results", required=True)
    parser.add_argument("--simul", type=int, help="number of iterations, default=0 (no siml)", default=0)
    parser.add_argument("--plot", help="default = no plot", action="store_const", const="Y", default="N")
    parser.add_argument("--xmin", type=float, help="min limit for x axis in density curves", default=0)
    parser.add_argument("--xmax", type=float, help="max limit for x axis in density curves", default=0.5)
    parser.add_argument("--skip_bed_dup", help="skip duplicated bed line (same seqid+start+end), default=Y", choices=["Y", "N"], default="Y")

    args = parser.parse_args()
    main() 
