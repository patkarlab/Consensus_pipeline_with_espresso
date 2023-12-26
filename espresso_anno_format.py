import sys
import pandas as pd

args = sys.argv
filename = args[1]
sample_name = args[2]
directory = args[3]
df = pd.read_csv(filename)

##GT:AD:VAF:Qual1:Qual2:MapQual1:MapQual2:Reads1Plus:Reads1Minus:Reads2Plus:Reads2Minus:FlankingSeqGroup:model:model_Pvalue:corrected_pvalue

x = df['Otherinfo1']

data = dict()

data.setdefault('REF_COUNT', [])
data.setdefault('ALT_COUNT', [])
data.setdefault('VAF', [])
data.setdefault('QUAL1', [])
data.setdefault('QUAL2', [])
data.setdefault('MAPQUAL1', [])
data.setdefault('MAPQUAL2', [])
data.setdefault('READ1PLUS', [])
data.setdefault('READ1MINUS', [])
data.setdefault('READ2PLUS', [])
data.setdefault('READ2MINUS', [])
data.setdefault('FLANKINGSEQ_GROUP', [])
data.setdefault('MODEL', [])
data.setdefault('PVALUE_ADJ', [])


for row in x:
    rowitems=row.split('/')
    formatval=rowitems[-1].split(':')
    RD_AD=formatval[1].split(',')

    readdepth=RD_AD[0] #RD
    altdepth=RD_AD[1] #AD
    vaf_value=float(formatval[2])
    #float(vaf_value)
    #print(vaf_value)
    vaf = vaf_value * 100
    qual1=formatval[3]
    #print(vaf)

    data['REF_COUNT'].append(readdepth)
    data['ALT_COUNT'].append(altdepth)
    data['VAF'].append(vaf)
    data['QUAL1'].append(formatval[3])
    data['QUAL2'].append(formatval[4])
    data['MAPQUAL1'].append(formatval[5])
    data['MAPQUAL2'].append(formatval[6])
    data['READ1PLUS'].append(formatval[7])
    data['READ1MINUS'].append(formatval[8])
    data['READ2PLUS'].append(formatval[9])
    data['READ2MINUS'].append(formatval[10])
    data['FLANKINGSEQ_GROUP'].append(formatval[11])
    data['MODEL'].append(formatval[12])
    data['PVALUE_ADJ'].append(formatval[13])

df1=df.iloc[:,:5]
df2=pd.DataFrame(data, columns=data.keys())
df3=df.iloc[:,5:]

horizontal_stack = pd.concat([df1, df2, df3], axis=1)
horizontal_stack['cosmic84']=horizontal_stack['cosmic84'].astype(str).str.replace(',' , ';')
horizontal_stack['AAChange.refGene']=horizontal_stack['AAChange.refGene'].astype(str).str.replace(',' , ';')
horizontal_stack.replace(to_replace='.', value='-1', inplace=True)
horizontal_stack.rename(columns = {'Func.refGene':'Variant Site', 'ExonicFunc.refGene':'Variant Function'}, inplace = True)
outfile = directory + sample_name + "_espresso" + ".csv"
horizontal_stack.to_csv(f'{outfile}', index=False)
