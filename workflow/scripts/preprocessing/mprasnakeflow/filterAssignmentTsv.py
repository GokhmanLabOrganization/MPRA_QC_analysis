#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: martin.kircher@bihealth.de
:Date: *17.09.2013
"""

import sys
from collections import defaultdict

import click

## Assumes a barcode sorted file with three columns
## 1: barcode
## 2: assigned insert
## 3: quality metrics of assignment (i.e. alignment)


def mfreq(datadict: defaultdict):
    mvalue = max(datadict.values())
    for key, value in datadict.items():
        if mvalue == value:
            return key
    return "NA"


def filter_and_report_assignments(
    fraction: float,
    minimum: int,
    ambiguous: bool,
    other: bool,
    cbarcode: str | None,
    cassignments: defaultdict,
    cquality: defaultdict,
):
    total = float(sum(cassignments.values()))
    if total == 0:
        return
    reported = False
    maxcount = 0
    insert = 0
    for insert, count in cassignments.items():
        if count > maxcount:
            maxcount = count
        else:
            continue
        if (total >= minimum) and (count >= minimum) and (count / total >= fraction):
            reported = True
            if (insert != "other") or (other):
                sys.stdout.write("%s\t%s\t%s\t%d/%d\n" % (cbarcode, insert, mfreq(cquality[insert]), count, total))
    if ambiguous and (not reported):
        sys.stdout.write("%s\t%s\t%s\t%d/%d\n" % (cbarcode, "ambiguous", mfreq(cquality[insert]), maxcount, total))


@click.command()
@click.option("-m", "--minimum", default=3, type=int, help="Minimum support for an assignment to get reported.")
@click.option("-f", "--fraction", default=0.75, type=float, help="Fraction of reads required to support assignment.")
@click.option("-o", "--other", is_flag=True, default=False, help="Also report barcodes for unknown/other insert sequence.")
@click.option("-a", "--ambiguous", is_flag=True, default=False, help="Also report barcodes for ambiguous insert sequence.")
def main(minimum: int, fraction: float, other: bool, ambiguous: bool):
    if fraction <= 0.5:
        raise click.BadParameter("Fraction has to be above 0.5.", param_hint="'--fraction'")

    cbarcode = None
    cassignments = defaultdict(int)
    cquality = defaultdict(lambda: defaultdict(int))

    for line in sys.stdin:
        fields = line.rstrip().split("\t")
        if len(fields) == 3:
            tag = fields[0]
            if tag != cbarcode:
                filter_and_report_assignments(fraction, minimum, ambiguous, other, cbarcode, cassignments, cquality)
                if (cbarcode is not None) and (tag < cbarcode):
                    raise click.ClickException("File is not sorted by barcode!")
                cbarcode = tag
                cassignments = defaultdict(int)
                cquality = defaultdict(lambda: defaultdict(int))
            cassignments[fields[1]] += 1
            cquality[fields[1]][fields[2]] += 1

    # final reporting
    filter_and_report_assignments(fraction, minimum, ambiguous, other, cbarcode, cassignments, cquality)


if __name__ == "__main__":
    main()
