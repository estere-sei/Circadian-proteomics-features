import re, sys, os
import numpy as np

from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rcParams['font.sans-serif'] = "Arial" 

PROG_NAME = 'extract_kinase_motifs'
DESCRIPTION = 'Identify protein kinase motifs, as collected at the PHOSIDA database, in FASTA files'
DEFAULT_SUFFIX = '_kinase_motifs'

KINASE_MOTIFS = {'PKA'     :[('R.[ST]', 2),
 		             ('R[RK].[ST]', -1),
                             ('KR..[ST]', -1)],
                 'CK1'     :[('S..[ST]', -1),
                             ('[ST]...S', -1)],
                 'CK2'     :[('[ST]..E', 0)],
                 'GSK3'    :[('S...S', 0)],
                 'CDK2'    :[('[ST]P.[KR]', 0)],
                 'CAMK2'   :[('R..[ST]', -1),
                             ('R..[ST]V', -2)],
                 'ERK/MAPK':[('P.[ST]P', 3),
                             ('V.[ST]P', 3),
                             ('PE[ST]P', 3)],
                 'PKB/AKT' :[('R[RST].[ST].[ST]', -3),
                             ('R.R..[ST]', -1)],
                 'PKC'     :[('R..[ST].R', -3)],
                 'PKD'     :[('[LVI].[RK]..[ST]', -1)],
                 'LCK'     :[('[IEV]Y[EG][EDPN][IVL]', 1)],
                 'ABL'     :[('[IVL]Y..[PF]', 1)],
                 'SRC'     :[('[ED]..Y..[DEAGST]', 4)],
                 'ALK'     :[('Y..[ILVM]', 0)],
                 'EGFR'    :[('[DPSAEN].Y[VLDEINP]', 2)],
                 'CDK1'    :[('[ST]P.[KR]', 0),
                             ('[ST]P[KR]', 0)],
                 'AURORA'  :[('[RK].[ST][ILV]', 2)],
                 'AURORA-A':[('[RKN]R.[ST][MLVI]', 3)],
                 'PLK'     :[('[DE].[ST][VILM].[DE]', 2)],
                 'PLK1'    :[('[ED].[ST][FLIYWVM]', 2)],
                 'NEK6'    :[('L..[ST]', -1)],
                 'CHK1/2'  :[('L.R..[ST]', -1)],
                 'CHK1'    :[('[MILV].[RK]..[ST]', -1)],
                 'PDK1'    :[('F..F[ST][FY]', -2)],
                 'NEK1'    :[('[FLM][RK][RK][ST]', -1)],
                }


def iter_fasta(stream_or_path):
  
  if isinstance(stream_or_path, str):
    stream = open(stream_or_path)
  else:
    stream = stream_or_path
    
  name = None
  seq = []
  
  for line in stream:
    line = line.strip()
    
    if not line:
      continue
    
    if line[0] == '>':
      if name:
        yield name, ''.join(seq)

      seq  = []
      name = line[1:]
    else:
      seq.append(line)

  if name:
    yield name, ''.join(seq)
         
                
def extract_kinase_motifs(fasta_paths, out_suffix=DEFAULT_SUFFIX, pdf_path=None, screen_gfx=False):
  
  motif_dict = {}
  for motif_name in KINASE_MOTIFS:
    motif_dict[motif_name] = []
    
    for site, offset in KINASE_MOTIFS[motif_name]:
      pattern = re.compile(site)
      motif_dict[motif_name].append(pattern)
  
  if not pdf_path:
    pdf_path = os.path.join(os.path.dirname(fasta_paths[0]), 'motif_counts.pdf')
  
  fig, ax = plt.subplots(figsize=(12, 5))
  
  motif_keys = sorted(motif_dict)
  n_keys = len(motif_keys)

  n_fasta = float(len(fasta_paths))
     
  from colorsys import hsv_to_rgb
  colors = [hsv_to_rgb(h, 1.0, 0.9) for h in np.arange(0.0, 0.8, 1.0/n_fasta)] 
  colors = ['#%02X%02X%02X' % (int(r*255), int(g*255), int(b*255)) for r,g,b in colors]
  bar_width = 0.8/n_fasta
  
  for f, fasta_path in enumerate(fasta_paths):
    fasta_name = os.path.splitext(os.path.basename(fasta_path))[0]
    out_path = '%s%s.tsv' % ( os.path.splitext(fasta_path)[0], out_suffix)
    
    with open(fasta_path) as file_obj, open(out_path, 'w') as out_file_obj:
      write = out_file_obj.write
      
      line = '\t'.join([''] + motif_keys) + '\n'
      write(line)
      bar_data = np.zeros(n_keys)
      
      out_data = []
      n_seqs = 0
      j = 0
      
      for name, seq in iter_fasta(file_obj):
        motifs = defaultdict(int)
        n_seqs += 1
        
        for motif_name in motif_keys:
          for pattern in motif_dict[motif_name]:
            match = pattern.search(seq)
            
            if match:
              motifs[motif_name] += 1
              break
     
        if motifs:
          j += 1
          
          line_data = [name] + ['%d' % motifs[x] for x in motif_keys]
          line = '\t'.join(line_data) + '\n'
          write(line)
          
          for k,  key in enumerate(motif_keys):
            bar_data[k] += motifs[key]
          
      msg = 'Processed {:,} sequences in {}: found {:,} kinase matches'
      print(msg.format(n_seqs, fasta_path, j))    
      
      msg = 'Written text file {}'
      print(msg.format(out_path))
      
      bar_data *= 100.0/n_seqs
      
      off = 0.1 + (0.8 * f)/n_fasta
      
      x_vals = np.arange(off, n_keys+off)
      
      ax.bar(x_vals, bar_data, width=bar_width, color=colors[f], label=fasta_name) 
 
 
  ax.set_xticklabels(motif_keys, fontsize=12, rotation=-90)
  ax.xaxis.set_ticks(np.arange(0.5, n_keys+0.5))
  ax.set_xlabel('Kinase site')
  ax.set_ylabel('% input sequences')
  
  plt.tight_layout()
  
  if screen_gfx:
    pdf = None
  else:
    pdf = PdfPages(pdf_path)
  
  ax.legend(fontsize=9)
  
  if pdf:
    pdf.savefig(dpi=100)
    plt.close()
    pdf.close()
    print('Written PDF file {}'.format(pdf_path))
  else:
    plt.show() 
    print('Done')   
    
           
if __name__ == '__main__':

  from argparse import ArgumentParser
 
  argv = sys.argv[1:]

  epilog = 'For further help email tjs23@cam.ac.uk'

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True ) 
                             
  arg_parse.add_argument(metavar='FASTA_FILES', nargs='+', dest='i',
                         help='Input one or more FASTA format files containing named protein sequences. Wildcards accepted.')
                         
  arg_parse.add_argument('-t', metavar='OUT_FILE_SUFFIX', default=DEFAULT_SUFFIX,
                         help='Optional suffix for naming tab-separated text files based in the input FASTA file. Default: %s' % DEFAULT_SUFFIX)

  arg_parse.add_argument('-o', '--out-pdf', metavar='PDF_FILE', default=None, dest="o",
                         help='Output PDF format file. If not specified, a default based on the input file name(s).')


  arg_parse.add_argument('-g', '--gfx', default=False, action='store_true', dest="g",
                         help='Display graphics on-screen using matplotlib and do not automatically save output.')

 
  args = vars(arg_parse.parse_args(argv))

  in_paths = args['i']        
  out_suffix = args['t']
  pdf_path = args['o']
  screen_gfx = args['g']
                         
  extract_kinase_motifs(in_paths, out_suffix, pdf_path, screen_gfx)

  # Bar chart of overall counts/kinase
