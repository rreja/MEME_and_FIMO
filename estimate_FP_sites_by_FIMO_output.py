
def process_file(infile,fimoFiles,options):
    input = open(infile,'rt')
    for line in input:
        if line.startswith("#"):
            continue
        if(len(line.split("\t")) == 9):
            countLines = countLines+1
            chrom,junk,junk,start,end,tag,strand,junk,attr = line.split("\t")
            fullfile[line] = float(tag)
    sortedfile = sorted(fullfile.iteritems(),key=operator.itemgetter(1),reverse=True)
    for k in sortedfile:
        bedfile = bedfile+k[0]
    for k,v in storedir.items():
        #outdir = os.path.join(os.path.dirname(fname),v)
        #outfile = os.path.join(outdir,os.path.splitext(os.path.basename(fname))[0]+"_"+v+".fasta")
        #outfile2 = os.path.join(outdir,os.path.splitext(os.path.basename(fname))[0]+"_"+v+".gff")
        retrieve_fasta(options.method,bedfile,k,outfile,countLines,maxocc,ref,options,outfile2)
        a = pybedtools.BedTool(bedfile, from_string=True).head(n=thresh,as_string=True)
        pybedtools.BedTool(a,from_string=True).slop(g=options.gfile,l=options.up_dist,r=options.down_dist).sequence(fi=ref,fo=outfile)
    










usage = '''
input_paths may be:
- a directory to run on all files in them
- a single file.


example usages:
python estimate_FP_sites_by_FIMO_output.py -i /usr/local/peak-pair.gff /usr/local/data/cwpair.gff -g sg07.txt -u 50 -d 50 --threshold="20 40 60" /usr/local/peak-pairs
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  

def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    #parser.add_option('-r', action='store', type='string', dest='peakpairFile',
    #                  help='Peak-pair mid point file in gff format.')
    parser.add_option('-f', action='store', type='string', dest='FIMOoutdir',
                      help='Dir cotaining FIMO output of various motifs from the given peak-pair file in gff format.')
    #parser.add_option('-r', action='store',  dest='reffile',
    #                  help='Reference genome in fasta format. Chromosme should match with input file')
    parser.add_option('-u', action='store', type='int', dest='up_dist',default=20,
                      help='No of bases to subtract from the start coordinate, Default=20')
    parser.add_option('-d', action='store',  type='int', dest='down_dist',default=20,
                      help='No of bases to add to the end coordinate, Default=20')
    parser.add_option('-g', action='store',  type='string', dest = 'gfile',
                      help='text file containing the size of each chromosome in the format: chr start end, ex: chr1 1 123456')
    
    (options, args) = parser.parse_args()
    
    
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
     
    # Create array to store all the fimo output files corresponding to different motifs.
    fimoFiles = []
    for name in os.listdir(options.FIMOoutdir):
            if name.endswith('.gff'):
                fimoFiles.append(os.path.join(options.FIMOoutdir, name)) 
        
    if not os.path.isdir(args[0]):
        process_file(args[0],fimoFiles,options)
        
    #else:
    #            
    #    if not os.path.exists(args[0]):
    #        parser.error('Path %s does not exist.' % args[0])
    #        
    #    for name in os.listdir(args[0]):
    #        if name.endswith('.gff'):
    #            fname = os.path.join(args[0], name)
    #            process_file(fname,options,storedir)
       

if __name__ == "__main__":
    run() 