#!/usr/bin/python

from math import sqrt

def spearmans(data_filename):
    
    ranks, N = collectRunningRanks(data_filename)
    
    scorr = calculateSpearmans(ranks, N)
    
    return scorr
### end spearmans


def pearsons(data_filename):
    
    """
    This program calculates the Pearsons moment correlation and
    Spearmans rank correlation for a set of columns in a lengthy file,
    meaning there are a lot of values in each column. The file is read
    line-by-line and intermediate statistics are updated that are finally
    used to calculate the correlations.

    The Pearson's correlation is calculated by:


    \rho_{X,Y}=
        \frac
        {
            \operatorname{E}[XY]-\operatorname{E}[X]\operatorname{E}[Y]
        }
        {
            \sqrt{\operatorname{E}[X^2]-\operatorname{E}[X]^2}
            ~
            \sqrt{\operatorname{E}[Y^2]- \operatorname{E}[Y]^2}
        }

    The expectations (moments) in the above equation are sums, and can
    be calculated cumulatively as the file is read line-by-line.

    The above can also be expressed as:

    r_{xy} = 
        \frac
        {
            N\sum x_iy_i-\sum x_i\sum y_i
        }
        {
            \sqrt
            {
                N\sum x_i^2-(\sum x_i)^2
            }
            ~
            \sqrt
            {
                N\sum y_i^2-(\sum y_i)^2
            }
        }

    The sample size, $N$ is the number of lines read,
    and does not need to be used until the very end when
    the correlation coefficient is calculated from the sums.
    """
    
    E_, E_XY, N = \
        collectRunningMoments(data_filename)
    
    pcorr = calculatePearsons(E_, E_XY, N)
    
    return pcorr
### end pearsons


def pileup_corr(data_filename, max_depth):
    
    """spearmans, pearsons = pileup_correlations(data_filename)"""
    
    # Collect the intermediate data
    ranks, E_, E_XY, N = \
        collectRunningMoments(data_filename)
        
    # Calculate the correlation coefficients
    pearsons = calculatePearsons(E_, E_XY, N)
    spearmans = calculateSpearmans(ranks, N)
    
    # Write the correlation matrices to files
    #
    # Another approach would be to pass in WHICH one you want
    # and send it stdout
    writeCorrelationFile(spearmans_correlation_file_name, spearmans)
    writeCorrelationFile(pearsons_correlation_file_name, pearsons)

    # return them in case called by other program
    return spearmans, pearsons
    
### end main
    

def collectRunningMoments(data_filename):

    """
    ranks, E_, E_XY, N =  collectRunningMoments(data_filename)
    
    Reads through the file line by line and collects the running sums
    and ranks which will be used to calculate the correlation coeeficients.
    
    Also keeps a running count of the number of values.
    
    """
    
    # get the number of sequences in the multipilup 
    # in order to initialize the correlation matrices
    num_columns = getNumColumns(data_filename)

    # initialize array of means
    E_ = zeros(num_columns)

    # Initialize the 2D arrays for correlation coeeficients and moment products
    pearsons = E_XY = zeros(num_columns, num_columns)
            
    # initialize the values from a single line of the multipileup file
    values = zeros(num_columns)
    
    # Read in the lines one by one and update the intermediate
    # stats. 
    N = 0
    with open(data_filename, 'r') as pileup_file:
        for line in pileup_file:
                        
            # 1) Split the line by whitespace.
            # 2) Skip the first three elements, then take the first
            #    element of every subsequent group of three.
            # See docstring for getNumColumns.
            
            # print("line " + line)
            values = [float(i) for i in line.split()]
            if len(values) != num_columns:
                break
            # Calculate moments. E_XX is on the diagonal of
            # E_XY.
            for i in range(num_columns):
                E_[i] += values[i]
                for j in range(num_columns):
                    E_XY[i][j] += values[i]*values[j]

            # increment number of values for normalization 
            N += 1
        ### end for
    ### end with
            
    return E_, E_XY, N
    
### end collectRunningMoments


def collectRunningRanks(data_filename):


    """
    ranks, N =  collectRunningRanks(data_filename)
d
    Reads through the file line by line and collects the ranks
    which will be used to calculate the correlation coeeficients.

    To calculate the rank we first need to count up the number
    of times each depth occurs. These counts will get sorted
    and processed to remove ties.

    Initialize a list containing a list of depth_counts for
    each column, as well as a max_depth and the number of
    values (N).

    The 1st dimension of depth_counts is the number
    of columns, and the 2nd dimension is the maximum depth.
    For instance, the number of times a depth of 29 occurs among
    the volues in column 5 is depth_counts[5][29].

    """

    # get the number of sequences in the multipilup 
    # in order to initialize the correlation matrices
    num_columns = getNumColumns(data_filename)
            
    # initialize the values from a single line of the multipileup file
    values = zeros(num_columns)

    # initialize depth_counts with empty lists
    depth_counts = [ [] for dummy in range(num_columns) ]
    max_depth = 0

    # pileup_file = open(data_filename, 'r')

    # Read in the lines one by one and update the intermediate
    # stats. 
    N = 0
    with open(data_filename, 'r') as data_file:
        for line in data_file:
            
            values = line.split()
            depth_counts = updateDepthCounts(values, depth_counts)
            
            N += 1
        ### end for
    ### end with
            
    return ranks, E_, E_XY, N

### end collectRunningMoments


def updatePearsonsMoments(values, E_, E_XY):

    """
    E_, X_XY = updatePearsonsMoments(values, E_, E_XY)

    Updates the running sums of the means and cross-products by 
    adding the values from the current line.
    """

    # note: len(values) should be equal to
    # num_columns, but that is not in scope and seems
    # awkward to pass in

    num_columns = len(values)

    for i in range(num_columns):
        E_[i] += values[i]
        for j in range(num_columns):
            E_XY[i][j] += values[i]*values[j]
        ### end j
    ### end i

    return E_, E_XY

### end updatePearsonsMoments


def calculatePearsons(E_, E_XY, N):
    """
    pearsons = calculatePearsons(E_, E_XY)

    E_      : 1xM sum of values
    E_XY    : MxM sum of products
    N       : Number of values
    """  

    num_columns = len(E_)
    cov = zeros(num_columns, num_columns)

    for i in range(num_columns):
        for j in range(num_columns):
            cov[i][j] = \
                (N*E_XY[i][j] - E_[i]*E_[j])/sqrt(N*E_XY[i][i] - E_[i]*E_[i])/sqrt(N*E_XY[j][j] - E_[j]*E_[j])
            #corr_coeff[i][j] = sqrt(cov[i][j])
        ### end j
    ### end i
    return cov

### end calculateSpearmans


def calculateSpearmans(ranks, N):
    """
    """

### end calculatePearsons


def updateDepthCounts(values, depth_counts):

    num_columns = len(values)
    # max_depth = max( [ max(depth_counts[i]) for i in range(num_columns) ] )

    # manage static max depth
    # first test if it exists
    try: test = updateDepthCounts.max_depth
    # if not, initialize to 0
    except AttributeError:  updateDepthCounts.max_depth = 0

    for i in range(num_columns):
        
        current_depth = values[i]
        
        if current_depth > updateDepthCounts.max_depth:
            # extend depth_counts[j] up to max_length for each j
            for i in range(num_columns): 
                depth_counts[i] += zeros(current_depth - updateDepthCounts.max_depth)
                
                # increment depth count
                depth_counts[i][current_depth] += 1
                
                # update max depth static variable
                updateDepthCounts.max_depth = current_depth

### end updateDepthCounts


def calculateRanks(values, ranks):

    """
    We need to calculate the ranks from the depth counts.
    Determine ranks from depth counts:

    1) sort counts
    2) 
    """


    ### end updateRanks


def calculateTiedRank(current_rank, depth):

    """

    sum 1...N = N*(N-1)/2
    sum N...M = M(M-1)/2 - (N-1)*(N-2)/2

    """

    N = current_rank + 1
    M = N + depth
    sum_of_raw_ranks = M*(M-1)/2 - (N-1)*(N-2)/2
    average_rank = sum_of_raw_ranks/depth
    return average_rank

### end calculateTiedRank


def zeros(N, M=0):
    """
    Z = zeros(N,M)

    if M is blank or 0, Z is a 1xN list

    if M is given, then Z is a NxM list
    """

    if M == 0:
        return [0 for n in range(N)]
    else:
        return [ [0 for m in range(M)] for n in range(N)]
    
### end zeros


def getNumColumns(data_filename):

    data_file = open(data_filename, 'r')
    line = data_file.readline()
    data_file.close()

    return len(line.split())

### end getNumColumns



def  writeCorrelationFile(correlation_file_name, corr):

    import sys

    N = len(corr)
    M = len(corr[0])
    
    # print(corr)

    if correlation_file_name != "":
        output_file = open(correlation_file_name, 'w')

    if N != M:
        print("The correlation matrix does not appear to be square. Something went wrong.")
        sys.exit()
    for i in range(N):
        row = " ".join([str(x) for x in corr[i]])
        print(row)
        if correlation_file_name != "":
            output_file.write(row + "\n")
    if correlation_file_name != "":
        output_file.close()
            

### end writeCorrelationFile



#~ for line in file of values at each position
#~ {
#~ value = parse out the depth;
#~ values[depth]++;
#~ }

#~ current_rank = 0;
#~ for (i = 0; i<values.length ; i++)
#~ {
#~ current_depth = values[i];
#~ if (current_depth == 0)
#~ {   
    #~ new_rank = current_rank;
#~ }
#~ else if (current_depth == 1)
#~ {
    #~ new_rank = floor(current_rank) + 1;
#~ }
#~ else
#~ {
    #~ // number of values with next highest rank = depth
    #~ // so their ordinal ranks are current_rank + 1 ...current_rank + 1 + depth
    #~ // the usual method is to assign them all a rank of their average 
    #~ // ordinal rank. In other words, the sum of their ranks/depth
    
    #~ new_rank = calculateTiedRank(current_rank, depth);
#~ }

#~ ranks[current_depth] = new_rank;
#~ }


if __name__ == "__main__":
    import sys

    args = sys.argv
    
    # print("args: " + ", ".join(args))

    if len(args) == 1:
        num_columns = len(values)
        print("ERROR: No arguments given.")
        print("Usage: corr_on_the_run corrtype input_file output_file")
        print("corrtype must be \"spearmans\" or \"pearsons\"")
        sys.exit()
        
    if len(args) == 2:
        print("ERROR: No input file given.")
        print("Usage: corr_on_the_run corrtype input_file output_file")
        print("corrtype must be \"spearmans\" or \"pearsons\"")
        sys.exit()
        
    if len(args) >= 3:
        corrtype = args[1]
        input_file = args[2]
        while True:
        
            # spearmans
            if corrtype == "spearmans":
                corr = spearmans(input_file)
                break
            
            # pearsons    
            if corrtype == "pearsons":
                corr = pearsons(input_file)
                break
                
            # help
            if corrtype == "help":
                print("Usage: corr_on_the_run corrtype input_file output_file")
                print("corrtype must be \"spearmans\" or \"pearsons\"")
                sys.exit()
                
            # corrtype error
            print("ERROR: corrtype must be \"spearmans\" or \"pearsons\"")
            print("Usage: corr_on_the_run corrtype input_file output_file")
            sys.exit()
            
        ##### end while
        
    ## end if

    if len(args) >= 4:
        output_file = args[3]
    else:
        output_file = ""
    
    writeCorrelationFile(output_file, corr)
