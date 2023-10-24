import sys
import os

help_flags = ['--help', '-h']
transcription_flags = ['--transcription', '-tc']
inverse_flags = ['--inverse', '-i']
translation_flags = ['--translation', '-tl']
verbose_flags = ['--verbose', '-v']
flags = transcription_flags + inverse_flags + translation_flags + verbose_flags

is_transcript = False
is_inverse = False
is_translate = False
is_verbose = False

extensions = ['.fasta', '.fna', '.ffn', '.faa', '.frn']

def print_help():
  print(""">>>>Transcriptor Script<<<<
usage: python bioScript.py <route> [ -h help | -v verbose | -tc transcription | -rtc retrotranscription | -tl translation ]
Options and arguments:
  -h,   --help \t\t\tShow this help message and exit.
  -v,   --verbose\t\tShow processes messages.
  -i,   --inverse\t\tSave the inverse chain of sequences (3'-5').
  -tc,  --transcription\t\tProcess chain only for transcription or retrotranscription (5'-3').
  -tl,  --translation\t\tProcess RNA chain only for translation.""")
  exit()

def print_error(message):
  print('\033[91m' + message + '\033[0m')
  print_help()
  
def print_warning(message):
  print('\033[93m' + message + '\033[0m')

def show(message):
  if is_verbose:
    print(message)

def verify_path(path):
  filename = os.path.splitext(os.path.basename(path))
  
  if filename[1] not in extensions:
    print_error("ERROR: file should be .fasta format")
  if not os.path.isfile(path):
    print_error("ERROR: please enter a valid file path")
  
  return filename 
  
args = sys.argv[1:]
n_flags = len(args)

route = ""
filename = []

if n_flags == 0:
  print_error("ERROR: There is no route specified")
elif n_flags == 1:
  if args[0] in help_flags: 
    print_help()
  if args[0] in flags:
    print_error("ERROR: There is no route specified")
    
  route = args[0]
  filename = verify_path(route)
  is_translate = is_transcript = True
 
else:
  for arg in args:
    if arg in flags:
      if arg in transcription_flags:
        is_transcript = True
      elif arg in translation_flags:
        is_translate = True
      elif arg in inverse_flags:
        is_inverse = True
      elif arg in verbose_flags:
        is_verbose = True
    else:
      route = arg
      filename = verify_path(route)
  
  if route == '':
    print_error("ERROR: There is no route specified")
  if is_verbose and not is_transcript and not is_inverse and not is_translate:
    is_translate = is_transcript = True
  

#---------- 1 TRANSCRIPTION ---------------#

transcription_3_5 = {
  'DNA' : {'A':'U', 'C':'G', 'T':'A', 'G':'C'},
  'RNA' : {'A':'T', 'C':'G', 'U':'A', 'G':'C'}
}
transcription_5_3 = {
  'DNA' : {'A':'A', 'C':'C', 'T':'U', 'G':'G'},
  'RNA' : {'A':'A', 'C':'C', 'U':'T', 'G':'G'}
}

def transcript_3_5(chain):
  show('Retrotranscription process for chain ' + chain['header'][1:-1])
  
  transcript = []
  for base in chain['chain']:
    transcript.append(transcription_3_5[chain['chain_type']][base])
  
  return transcript

def transcript_5_3(chain):
  show('Transcription process for chain ' + chain['header'][1:-1])
  
  transcript = []
  for base in chain['chain']:
    transcript.append(transcription_5_3[chain['chain_type']][base])
  
  return transcript

#---------- 2 TRADUCTION ------------------#

aminoacid_trans = { 'AUG':'M', #MET
                   'UAA':'-', 'UAG':'-', 'UGA':'-', #STOP
                   'UUA':'L', 'UUG':'L', 'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L', #Leu
                   'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S', 'AGU':'S', 'AGC':'S', #Ser
                   'CGU':'R','CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R', #Arg
                   'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', #Gly
                   'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', #pro
                   'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', #tHR
                   'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', #ala
                   'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',#val
                   'AUU':'I', 'AUC':'I', 'AUA':'I', #Ile
                   'UAU':'Y', 'UAC':'Y', #Tyr
                   'UUU':'F', 'UUC':'F', # Phe
                   'CAU':'H', 'CAC':'H', #His
                   'CAA':'Q', 'CAG':'Q', #Gln
                   'AAU':'N', 'AAC':'N', #Asn
                   'AAA':'K', 'AAG':'K', #Lys
                   'GAU':'D', 'GAC':'D', #Asp
                   'GAA':'E', 'GAG':'E', #Glu
                   'UGU':'C', 'UGC':'C', #Cys
                   'UGG':'W' #Trp
                  }

def translate(chain, frame):
  show('Transcription process for chain ' + chain['header'][1:-1])
  
  start = frame - 1
  aminoacids = []
  for i in range(3+start,len(chain['chain'])+1, 3):
    aminoacid = ''.join(chain['chain'][start: i])
    aminoacids.append(aminoacid_trans[aminoacid])
    start = i
  return aminoacids

#---------- 3 SAVE FILE -------------------#

def save_file(filename, chains):
  print('Saving chains in file', filename)
  
  chainfile = open(filename, 'w')
  for chain in chains.keys():
    chainfile.write(chains[chain]['header'])
    
    chain_len = chains[chain]['len']
    line = ''
    count = 0
    for base in chains[chain]['chain']:
      if count == chain_len:
        chainfile.write(line + '\n')
        line = ''
        count = 0
      
      line = line + base
      count = count + 1
    
    if count > 0:
      chainfile.write(line + '\n')
  chainfile.close()

#---------- Multi Fasta -------------------#
fastafile = open(route, "r")

chains = {}
n_chains = 0
is_valid = True

line = fastafile.readline()

#Validation process

valid_bases = ['A', 'T', 'U', 'G', 'C']

while len(line) != 0:  
  if line[0] == '>':
    n_chains = n_chains + 1
    show('processing chain ' + line[1:-1])
    chains[n_chains] = {'header': line,'chain_type' : '', 'chain': [], 'len' : 0}
    is_valid = True
    
  elif is_valid:
    n_bases = [*line][:-1]
        
    i = 0
    while is_valid and i < len(n_bases):
      base = n_bases[i]
      i = i + 1
      if base == 'T' and chains[n_chains]['chain_type'] == '':
        chains[n_chains]['chain_type'] = 'DNA'
        show('\t Identified as DNA...')
      elif base == 'T' and chains[n_chains]['chain_type'] != 'DNA':
        print_warning('WARNING: ' + chains[n_chains]['header'][1:-1] +' not a valid RNA chain. SKipping to the next chain...')
        is_valid = False
        n_chains = n_chains -1
      elif base == 'U' and chains[n_chains]['chain_type'] == '':
        chains[n_chains]['chain_type'] = 'RNA'
        show('\t Identified as RNA...')
      elif base == 'U' and chains[n_chains]['chain_type'] != 'RNA':
        print_warning('WARNING: ' + chains[n_chains]['header'][1:-1] + ' not a valid DNA chain. SKipping to the next chain...')
        is_valid = False
        n_chains = n_chains -1
      elif base not in valid_bases:
        print_warning('WARNING: '+ base + ' not a valid nitrogenous base. SKipping to the next chain...')
        is_valid = False
        n_chains = n_chains -1
        
    if is_valid:
      if chains[n_chains]['len'] == 0:
        chains[n_chains]['len']  = len(n_bases)
      chains[n_chains]['chain'].extend(n_bases) 
      
  line = fastafile.readline()

fastafile.close()

# Transcription process

retro_chains = {}
for i in chains.keys():
  chain = chains[i]
  if is_inverse:
    transcript = transcript_3_5(chain)
    retro_chains[i] = {'header': chain['header'] ,'len': chain['len'] , 'chain': transcript}
  
  chain['chain'] = transcript_5_3(chain)

if is_transcript and is_inverse:
  save_file('tc_' + filename[0] + "_5_3" + filename[1], retro_chains)
  
if is_transcript:
  save_file('tc_' + filename[0] + filename[1], chains)    

n_chains = len(chains.keys())
for i in range(1,n_chains):
  chain = chains[i]
  if chain['chain_type'] == 'RNA':
    del chains[i]

# --- frame 1 ---
chains_f1 = {}
for i in chains.keys():
  aminoacids = translate(chains[i], 1)
  chains_f1[i] = {'header': chains[i]['header'], 'len': chains[i]['len'], 'chain': aminoacids}

if is_translate:
  save_file('tl_frame1_' + filename[0] + filename[1], chains_f1)
  
# --- frame 2 ---
chains_f2 = {}
for i in chains.keys():
  aminoacids = translate(chains[i], 2)
  chains_f2[i] = {'header': chains[i]['header'], 'len': chains[i]['len'], 'chain': aminoacids}

if is_translate:
  save_file('tl_frame2_' + filename[0] + filename[1], chains_f2)

# --- frame 3 ---
chains_f3 = {}
for i in chains.keys():
  aminoacids = translate(chains[i], 3)
  chains_f3[i] = {'header': chains[i]['header'], 'len': chains[i]['len'], 'chain': aminoacids}

if is_translate:
  save_file('tl_frame3_' + filename[0] + filename[1], chains_f3)
  
exit()