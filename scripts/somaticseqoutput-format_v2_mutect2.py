#!/usr/bin/env python3

import pandas as pd
import sys

args = sys.argv

filename = args[1]
outfile = args[2]
vartools=['MuTect2 | ', 'VarScan2 | ', 'VarDict | ', 'LoFreq | ', 'Strelka | ', 'Freebayes | ', 'Platypus | ']

df = pd.read_csv(filename)
x = df['Otherinfo1']
discarded_column=df.columns.get_loc('Otherinfo2')
data = dict()
somatic_cols=['Chr','Start','End','Ref','Alt','FILTER','SOMATIC_FLAG','REF_COUNT','ALT_COUNT','VAF%','REF_F2R1','REF_F1R2','ALT_F2R1','ALT_F1R2','Func.refGene','Gene.refGene','ExonicFunc.refGene','AAChange.refGene','Gene_full_name.refGene','Function_description.refGene','Disease_description.refGene','cosmic84','PopFreqMax']
data.setdefault('FILTER', [])
data.setdefault('SOMATIC_FLAG', [])
#data.setdefault('VariantCaller_Count', [])
#data.setdefault('Variant_Callers', [])
data.setdefault('VAF%', [])
data.setdefault('REF_COUNT', [])
data.setdefault('ALT_COUNT', [])
data.setdefault('REF_F2R1', [])
data.setdefault('REF_F1R2', [])
data.setdefault('ALT_F2R1', [])
data.setdefault('ALT_F1R2', [])
for row in x:
    #print (row)	
    rowitems=row.split('\t')
    data['FILTER'].append(rowitems[9])
    #print (rowitems[9])
    info=rowitems[10].split(';')
    formatval=rowitems[-2].split(':')
    #print (rowitems[-2].split(':'))
    #print (formatval)
    readcounts=formatval[1]
    #print (formatval[1])
    if len(info)!=5 and info[0]!='SOMATIC':
        data['SOMATIC_FLAG'].append('NON SOMATIC')
        name_toolskey=0
    else:
        data['SOMATIC_FLAG'].append(info[0])
        name_toolskey=1
    num_toolskey=name_toolskey+1    
    #tools_binary=info[name_toolskey].split('=')[1].split(',')
    #tools_name=''
    #for key in range(len(tools_binary)):
    #    if tools_binary[key]=='1':
    #        tools_name+=vartools[key]
    #data['Variant_Callers'].append(tools_name)        
    #data['VariantCaller_Count'].append(info[num_toolskey].split('=')[1])
    #vaf="{:.2%}".format(float(info[-1].split('=')[1]))
    #vaf="{:.2}".format(float(formatval[2])*100) 
    #data['VAF%'].append(vaf)
    data['REF_COUNT'].append(formatval[1].split(',')[0])
    data['ALT_COUNT'].append(formatval[1].split(',')[1])
    vaf_calc=(float(formatval[1].split(',')[1]) /( float(formatval[1].split(',')[0])  + float(formatval[1].split(',')[1]) )) * 100
    data['VAF%'].append(round (vaf_calc, 2 ))
    data['REF_F2R1'].append(formatval[-1])
    data['REF_F1R2'].append(formatval[-2])
    data['ALT_F2R1'].append(formatval[4])
    data['ALT_F1R2'].append(formatval[3])
            
df1=df.iloc[:,:5]
#print (data.keys())
df2=pd.DataFrame(data, columns=data.keys())
df3=df.iloc[:,5:discarded_column]

horizontal_stack = pd.concat([df1, df2, df3], axis=1)
horizontal_stack['cosmic84']=horizontal_stack['cosmic84'].str.replace(',' , ';')
horizontal_stack['AAChange.refGene']=horizontal_stack['AAChange.refGene'].str.replace(',' , ';')
horizontal_stack.replace(to_replace='.', value='-1', inplace=True)
horizontal_stack=horizontal_stack.reindex(columns = somatic_cols)
horizontal_stack.to_csv(outfile, index=False)
