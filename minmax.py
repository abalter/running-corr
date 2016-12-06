#!/usr/bin/env python
import sys

def findMinMax(filename, col, delimiter):
    
    min = sys.maxsize
    max = -sys.maxsize

    with open(filename, 'r') as f:
        for line in f:
            try:
                x = int(line.split(delimiter)[col - 1])
            except IndexError:
                print("ERROR: There are less than {} columns.".format(col))
                sys.exit()               

            min = x if x < min else min
            max = x if x > max else max
    
    return min, max

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Find max and min of a specific column',
        epilog="""A future enhancement will be to handle multiple columns."""
        )
    parser.add_argument('--file', '-f', 
        required=True, 
        type=str, 
        help="The input file"
        )
    parser.add_argument('--column', '-c', 
        required=True, 
        type=int, 
        help="The column number, 1 indexed like `cut` and other shell tools.",
        )
    parser.add_argument('--delimiter', '-d',
        required=False, 
        type=str, 
        help="Column separator. Defaults to tab.",
        default="\t"
        )
    
    args = parser.parse_args()

    min, max = findMinMax(args.file, args.column, args.delimiter)

    print("{}\t{}".format(min, max))

