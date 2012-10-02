import sys, os, pybedtools, operator
from optparse import OptionParser , IndentedHelpFormatter
import matplotlib.pyplot as plt



def process_file(infile,outfile1,outfile2,options):
    
    if options.strandedness == 0:
        gff_file = pybedtools.BedTool(infile).slop(g=options.gfile,l=options.up_dist,r=options.down_dist)
    else:
        gff_file = pybedtools.BedTool(infile).slop(g=options.gfile,l=options.up_dist,r=options.down_dist,s=True)
    
    blacklisted_gff = pybedtools.BedTool(options.BlacklistFile)
    gff_file.intersect(blacklisted_gff,v=True).saveas(outfile1)
    gff_file.intersect(blacklisted_gff,u=True).saveas(outfile2)


    

usage = '''
input_paths may be:
- a directory to run on all files in them
- a single file.


example usages:
python Remove_blacklisted_regions.py -b <black listed file in gff>  -g <Chr length file> -u 50 -d 50  <dir to your files>
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  

def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    #parser.add_option('-r', action='store', type='string', dest='peakpairFile',
    #                  help='Peak-pair mid point file in gff format.')
    parser.add_option('-b', action='store', type='string', dest='BlacklistFile',
                      help='File containing black listed regions in gff format.')
    parser.add_option('-u', action='store', type='int', dest='up_dist',default=0,
                      help='No of bases to subtract from the start coordinate, Default=0')
    parser.add_option('-d', action='store',  type='int', dest='down_dist',default=0,
                      help='No of bases to add to the end coordinate, Default=20')
    parser.add_option('-s', action='store',  type='int',dest='strandedness',default=0, 
                      help='If your file contains strand information, Default no strandedness, s=0, change to s=1 for strandedness')
    parser.add_option('-g', action='store',  type='string', dest = 'gfile',
                      help='text file containing the size of each chromosome in the format: chr start end, ex: chr1 1 123456.Required Parameter.')
    
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
     
    if not os.path.isdir(args[0]):
        # Create dir names for noS and onlyS folders
        outdir = os.path.join(os.path.dirname(args[0]),"_blacklist_removed")
    else:
        outdir = os.path.join(args[0],"_blacklist_removed")
      
    if not os.path.exists(outdir):
            os.mkdir(outdir)
        
    if not os.path.isdir(args[0]):
        outfile1 = os.path.join(outdir,os.path.splitext(os.path.basename(args[0]))[0]+"noB.gff")
        outfile2 = os.path.join(outdir,os.path.splitext(os.path.basename(args[0]))[0]+"onlyB.gff")
        process_file(args[0],outfile1,outfile2,options)
        
    else:
                
        if not os.path.exists(args[0]):
            parser.error('Path %s does not exist.' % args[0])
            
        for name in os.listdir(args[0]):
            if name.endswith('.gff'):
                fname = os.path.join(args[0], name)
                outfile1 = os.path.join(outdir,os.path.splitext(name)[0]+"noB.gff")
                outfile2 = os.path.join(outdir,os.path.splitext(name)[0]+"onlyB.gff")
                process_file(fname,outfile1,outfile2,options)
       

if __name__ == "__main__":
    run() 