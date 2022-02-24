import json
import urllib2
import sys
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--seqid_file", help="Input file containing the list of seqids")
parser.add_argument("--outdir", help="Location of the output directory")
args = parser.parse_args()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

fin = open(args.seqid_file)

for l in fin:
	seqid = l.strip()
	#data = 'seqids=["P04637"]'
	data = 'seqids=["'+seqid+'"]'
	request = urllib2.Request('http://d2p2.pro/api/seqid', data)
	response = json.loads(urllib2.urlopen(request).read())

	try:
		fout = open(args.outdir+"/"+seqid+'_d2p2_result.out', 'w')
		for record in response[seqid][0][2]['disorder']['consranges']:
			#print "\t".join(record)
			fout.write("\t".join(record)+"\n")
		fout.close()
	except:
		print "Sequence id :"+seqid+" has not been annotated in the database\n"
fin.close()
