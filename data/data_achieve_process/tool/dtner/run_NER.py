# -*- coding: UTF-8 -*-    
#!/usr/bin/env python
# label the possible name of databases and tools
import os 
import sys
import time
import dtNER_prepare
import glob
import json
import shutil

def run_NER(source_text):
	op_path = os.getcwd()
	tmp_stemp = "tmp_"+str(int(time.time()))
	tmp_path = os.path.join(op_path,tmp_stemp)
	os.mkdir(tmp_path)
	tmp_NER_input = os.path.join(tmp_path,"deploy.txt")
	print(sys.argv[0])
	main_ner_src = os.path.join(op_path,sys.argv[0].replace("run_NER.py",""))
	print(main_ner_src)
	main_ner_src = os.path.join(main_ner_src,"src")
	re_tag = dtNER_prepare.NER_prepare(source_text)
	with open(tmp_NER_input, 'w') as out_handle:
		for idx, term in enumerate(re_tag['word']):
			out_handle.write(re_tag['word'][idx]+" "+re_tag['pos'][idx]+" "+re_tag['chunk'][idx]+" "+re_tag['mor'][idx]+" "+re_tag['tag'][idx]+"\n")
	out_handle.close()
	os.chdir(main_ner_src)
	outside_order = 'python main.py --train_model=False --use_pretrained_model=True --dataset_text_folder=%s --output_folder=%s --pretrained_model_folder=../trained_models/bmisner' %(tmp_path,tmp_path)
	os.system(outside_order)
	result_path = glob.glob(os.path.join(tmp_path,"tmp*"))[0]
	result_file = os.path.join(result_path,'000_deploy.txt')
	re_tag = []
	with open(result_file,'r') as handle:
		for line in handle:
			a_line = line.strip()
			if a_line == '': continue
			a_line = a_line.split(" ")
			if(a_line[8]=="I-BME" or a_line[8]=="B-BME"): # new
				re_tag.append({"term":a_line[0],"type":a_line[8]}) #new
#			re_tag.append({"term":a_line[0],"F1":a_line[4],"F2":a_line[5], "F3":a_line[6], "F4": a_line[7], "Result": a_line[8] })
		shutil.rmtree(tmp_path, ignore_errors=True)
		# return(json.dumps(re_tag)) #new 
		return(re_tag) #new 


if __name__ == "__main__":
	kk=run_NER(sys.argv[1])
	print(kk)
	print(sys.argv[2])
	with open(sys.argv[2],"w") as fileout:
		for i in kk:
			fileout.write(i["term"]+"\t"+i["type"]+"\n")
