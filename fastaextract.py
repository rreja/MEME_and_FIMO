import sys, os
from scipy.stats import stats as st
from optparse import OptionParser , IndentedHelpFormatter
import pybedtools
import operator

def process_file(fname,options,storedir):
    input = open(fname,'rt')
    fullfile = {}
    bedfile = ""
    countLines = 0
    ref = pybedtools.BedTool(options.reffile)
    for line in input:
        if line.startswith("#"):
            continue
        if(len(line.split("\t")) == 9):
            countLines = countLines+1
            chrom,junk,junk,start,end,tag,strand,junk,attr = line.split("\t")
            fullfile[line] = float(tag)
    sortedfile = sorted(fullfile.iteritems(),key=operator.itemgetter(1),reverse=True)
    maxocc = float(sortedfile[0][1])
    for k in sortedfile:
        bedfile = bedfile+k[0]
    for k,v in storedir.items():
        outdir = os.path.join(os.path.dirname(fname),v)
        outfile = os.path.join(outdir,os.path.splitext(os.path.basename(fname))[0]+"_"+v+".fasta")
        retrieve_fasta(options.method,bedfile,k,outfile,countLines,maxocc,ref,options)
        
    
def retrieve_fasta(mthd,bedfile,thresh,outfile,totalLines,maxocc,ref,options):
    if(mthd == "count"):
        a = pybedtools.BedTool(bedfile, from_string=True).head(n=thresh,as_string=True)
        pybedtools.BedTool(a,from_string=True).slop(g=options.gfile,l=options.up_dist,r=options.down_dist).sequence(fi=ref,fo=outfile)
    elif(mthd == "pmo"):
        if thresh > 100:
            print "Your threshold for method=pmo should not be more than 100%"
        cutoff = float(thresh)/100*(maxocc)
        a = pybedtools.BedTool(bedfile, from_string=True)
        a.filter(lambda m:float(m.score) > cutoff).slop(g=options.gfile,l=options.up_dist,r=options.down_dist).sequence(fi=ref,fo=outfile)
    elif(mthd == "pct"):
        if thresh > 100:
            print "Your threshold for method=pct should not be more than 100%"
        cutoff = int((float(thresh)/100)*totalLines)
        a = pybedtools.BedTool(bedfile, from_string=True).head(n=cutoff,as_string=True)
        pybedtools.BedTool(a,from_string=True).slop(g=options.gfile,l=options.up_dist,r=options.down_dist).sequence(fi=ref,fo=outfile)
        

usage = '''
input_paths may be:
- a directory to run on all files in them
- a single file.


example usages:
python fastaextract -r /usr/local/ref.fasta /usr/local/data/cwpair.gff -g sg07.txt -u 50 -d 50 --threshold="20 40 60" /usr/local/peak-pairs
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  

def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('--thresholds', action='store', type='string', dest='th',default="10 50 100",
                      help='Types thresholds in "" seperated by single space, deafult "10 50 100"')
    parser.add_option('-m', action='store', type='string', dest='method',default="pct",
                      help='pmo=>perentage of max occupancy,cnt=>count,pct=> percent of all sorted entries,Default=pct')
    parser.add_option('-r', action='store',  dest='reffile',
                      help='Reference genome in fasta format. Chromosme should match with input file')
    parser.add_option('-u', action='store', type='int', dest='up_dist',default=0,
                      help='No of bases to subtract from the start coordinate, Default=0')
    parser.add_option('-d', action='store',  type='int', dest='down_dist',default=0,
                      help='No of bases to add to the end coordinate, Default=0')
    parser.add_option('-g', action='store',  type='string', dest = 'gfile',
                      help='text file containing the size of each chromosome in the format: chr start end, ex: chr1 1 123456')
    
    (options, args) = parser.parse_args()
    thresholds = options.th.split(" ")
    dirname = {}
    storedir = {}
    
    for k in thresholds:
        dirname["outdir_"+k] = "top"+k   # This stores the informtion of which filepointer corresponds to which directory.
        storedir[int(k)] =  "top"+k  # this stores which threshold corresponds to which file pointer.
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        
    try:
        from scipy.stats import stats
        import pybedtools
    except ImportError:
        print "You need to install Scipy(http://www.scipy.org/Download) before you can run this script."
        print "If you have Scipy installed, then you should check if you have pybedtools installed"
        sys.exit(1)
     
    #REF = open(options.reffile,'rt')    
    # Create direcotry that will store 1%, 10% etc etc information here for both cases.
    for k,v in dirname.items():
        k = os.path.join(os.path.dirname(args[0]),v)
        if not os.path.exists(k) :
            os.mkdir(k)
        
    if not os.path.isdir(args[0]):
        process_file(args[0],options,storedir)
        
    else:
                
        if not os.path.exists(args[0]):
            parser.error('Path %s does not exist.' % args[0])
            
        for name in os.listdir(args[0]):
            if name.endswith('.gff'):
                fname = os.path.join(args[0], name)
                process_file(fname,options,storedir)
       

if __name__ == "__main__":
    run() 