#i = raw_input("Put job number : ")

lines = []
total = 0
fake = 0
see = []

for i in range(50):
  with open("/data4/Users/jihkim/SKFlatRunlog/2019_10_16_205309__168954__ChargeFlipValidation__Year2016__CFrate__TAMSA1/DYJets/job_{0}.log".format(i),'r') as f:
    lines_tmp = f.readlines()
  
  for line in lines_tmp:
    if '!!' in line:
      total+=1
    if 'electrons :' in line or 'charge flipped electron pT' in line or 'charge flipped electron Eta, Phi' in line:
      lines.append(line[:-1])
    if 'Matched gen index' in line:
      try:
        gen_idx = int(line[-4:-1])
      except:
        gen_idx = int(line[-3:-1])
    try:
      int(line.split('\t')[0])
    except:
      pass
    else:
      if gen_idx == int(line.split('\t')[0]):
        lines.append('Matched gen PID : '+line.split('\t')[1])
        if abs(int(line.split('\t')[1])) != 11 :
          fake+=1
          lines.append('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
          see.append(i)
          see.append(lines[-3])
    

for j in range(len(lines)):
  print lines[j]

print 'CF events :', total
print 'fake electrons :', fake

if fake != 0:
  print 'See :'
  print see
 

 

