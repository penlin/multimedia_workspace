#coding:utf-8
import sys


SNR = 5

filename = str(sys.argv[-1])

fw = open('backup\\'+filename,'w')
fr = open(filename,'r')
data = [line.replace('\n','') for line in fr if line[0]=='0']
print 'Len:',len(data)
height = len(data)/SNR

psnr = [[' ' for j in range(SNR)] for i in range(height)]

ii = 0
jj = 0
for line in data:
    psnr[ii][jj] = line.split(', ')[-1]+', '
    ii=ii+1
    if ii==height :
        ii = 0
        jj = jj+1

for i in range(height):
    for j in range(SNR):
        fw.write(psnr[i][j])
    fw.write('\n')
		
fr.close()
fw.close()