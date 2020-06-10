# -*- coding: UTF-8 -*-    
#!/usr/bin/env python
# label the possible name of databases and tools
import os 
import sys
import time
import glob
import json
import shutil
import re
import nltk
import nltk.data
from nltk.corpus import conll2000

# for split text to sentence and sentence to words

def split_to_sentence(paragraph):
    tokenizer = nltk.data.load("tokenizers/punkt/english.pickle")
    sentences = tokenizer.tokenize(paragraph)
    return sentences

def split_to_word(sentence):
    pattern = r"""(?x)
                   (?:[A-Z]\.)+                  # abbreviations, e.g. U.S.A.
                   |(?:\d+[a-zA-Z]+)+
                   |\d+(?:[-._/]\d+)*%?           # numbers, incl. currency and percentages 
                   |\w+(?:[-'_]\w+)*             # words with optional internal hyphens
                   |\.\.\.                       # ellipsis 
                   |(?:[.,;/\|&!\[\]{}<>@#$^*+="'?~():-_`])        # special characters with meanings
                   """
    words = nltk.regexp_tokenize(sentence, pattern)
    return(words)
    
####### for morpholgy and por

def is_punc(a_word):
    all_punc = ".,;/\|&![]{}<>@#$%^*+=\"'?~():-_`"
    for letter in a_word:
        if letter not in all_punc:
            return False
    return True

def is_UL(a_word):
    if a_word[0].islower():
        return False
    else:
        for letter in a_word[1:]:
            if letter.isupper():
                return False
        return True

def is_AP(a_word):
    for letter in a_word:
        if (not letter.isalpha()) and (not is_punc(letter)):
            return False
    return True

def is_NPN(a_word):
    point_count = 0
    for letter in a_word:
        if letter == '.':
            point_count += 1
        elif not letter.isdigit():
            return False
    if point_count == 1:
        return True
    else:
        return False

def is_NP(a_word):
    for letter in a_word:
        if (not letter.isdigit()) and (not is_punc(letter)):
            return False
    return True

def have_upper_letter(a_word):
    for letter in a_word:
        if letter.isupper():
            return True
    return False


######################### for chunking

class UnigramChunker(nltk.ChunkParserI):
    """
    一元分块器，
    该分块器可以从训练句子集中找出每个词性标注最有可能的分块标记，
    然后使用这些信息进行分块
    """

    def __init__(self, train_sents):
        """
        构造函数
        :param train_sents: Tree对象列表
        """
        train_data = []
        for sent in train_sents:
            # 将Tree对象转换为IOB标记列表[(word, tag, IOB-tag), ...]
            conlltags = nltk.chunk.tree2conlltags(sent)

            # 找出每个词性标注对应的IOB标记
            ti_list = [(t, i) for w, t, i in conlltags]
            train_data.append(ti_list)

        # 使用一元标注器进行训练
        self.__tagger = nltk.UnigramTagger(train_data)

    def parse(self, tokens):
        """
        对句子进行分块
        :param tokens: 标注词性的单词列表
        :return: Tree对象
        """
        # 取出词性标注
        tags = [tag for (word, tag) in tokens]
        # 对词性标注进行分块标记
        ti_list = self.__tagger.tag(tags)
        # 取出IOB标记
        iob_tags = [iob_tag for (tag, iob_tag) in ti_list]
        # 组合成conll标记
        conlltags = [(word, pos, iob_tag) for ((word, pos), iob_tag) in zip(tokens, iob_tags)]

        return conlltags




def split_and_tag(origin_text):
    re_tag = {}
    re_tag["word"]=[]
    re_tag["tag"] = []
    a_line = origin_text.strip()
    sentences = split_to_sentence(a_line)
    sentences_words = [split_to_word(a_sentence) for a_sentence in sentences]
    for a_sentence in sentences_words:
        for a_word in a_sentence:
            res_line = a_word + '\t' + 'O' + '\n'
            re_tag["word"].append(a_word)
            re_tag["tag"].append("O")
    return(re_tag)
    


def get_pos_and_re_morphology(re_tag):
    re_tag["pos"] = []
    re_tag["mor"] = []
    words_length = len(re_tag['word'])
    words_pos = nltk.pos_tag(re_tag['word'])
    for a_word_pos in words_pos:
        re_tag["pos"].append(a_word_pos[1])
    for a_word in re_tag['word']:
        if a_word.isalpha():
            if a_word.isupper():
                if len(a_word) == 1:
                    re_tag["mor"].append('U')
                else:
                    re_tag["mor"].append("UU")
            elif a_word.islower():
                re_tag["mor"].append("LL")
            else:
                if is_UL(a_word):
                    re_tag["mor"].append("UL")
                else:
                    if a_word[0].isupper():
                        re_tag["mor"].append("ULU")
                    else:
                        re_tag["mor"].append("LUL")
        elif a_word.isdigit():
            if len(a_word) == 1:
                re_tag["mor"].append('N')
            else:
                re_tag["mor"].append("NN")
        elif a_word == '(':
            re_tag["mor"].append("LB")
        elif a_word == ')':
            re_tag["mor"].append("RB")
        elif a_word == ':':
            re_tag["mor"].append("CO")
        elif is_punc(a_word):
            re_tag["mor"].append("OP")
        elif a_word.isalnum():
            if have_upper_letter(a_word):
                re_tag["mor"].append("UAN")
            else:
                re_tag["mor"].append("LAN")
        elif is_AP(a_word):
            if have_upper_letter(a_word):
                re_tag["mor"].append("UAP")
            else:
                re_tag["mor"].append("LAP")
        elif is_NPN(a_word):
            re_tag["mor"].append("NPN")
        elif is_NP(a_word):
            re_tag["mor"].append("NP")
        else:
            if have_upper_letter(a_word):
                re_tag["mor"].append("UANP")
            else:
                re_tag["mor"].append("LANP")
    return(re_tag)




def nltk_chunking(re_tag):
    train_sents = conll2000.chunked_sents("train.txt")
    unigram_chunker = UnigramChunker(train_sents)
    re_tag["chunk"]=[]
    a_sentence = []
    for idx, term in enumerate(re_tag['word']):
      a_sentence.append((re_tag['word'][idx],re_tag['pos'][idx]))
    sentence_length = len(a_sentence)
    if sentence_length > 0:
        chunking_sentence = unigram_chunker.parse(a_sentence)
        for i in range(sentence_length):
            chunking_list = list(chunking_sentence[i])
            if chunking_list[-1] == None:
                chunking_list[-1] = 'O'
            re_tag['chunk'].append(chunking_list[-1])
    return(re_tag)


def NER_prepare(origin_text):
    re_tag = split_and_tag(origin_text)
    re_tag = get_pos_and_re_morphology(re_tag)
    re_tag = nltk_chunking(re_tag)
    return(re_tag)

def run_NER(source_text):
	op_path = os.getcwd()
	print(op_path)
	tmp_stemp = "tmp_"+str(int(time.time()))
	tmp_path = os.path.join(op_path,tmp_stemp)
	os.mkdir(tmp_path)
	tmp_NER_input = os.path.join(tmp_path,"deploy.txt")
	print(sys.argv[0])
	main_ner_src = os.path.join(op_path,"tool/dtner/")
	print(main_ner_src)
	main_ner_src = os.path.join(main_ner_src,"src")
	re_tag = NER_prepare(source_text)
	with open(tmp_NER_input, 'w') as out_handle:
		for idx, term in enumerate(re_tag['word']):
			out_handle.write(re_tag['word'][idx]+" "+re_tag['pos'][idx]+" "+re_tag['chunk'][idx]+" "+re_tag['mor'][idx]+" "+re_tag['tag'][idx]+"\n")
	out_handle.close()
	os.chdir(main_ner_src)
	outside_order = 'python main.py --train_model=False --use_pretrained_model=True --dataset_text_folder=%s --output_folder=%s --pretrained_model_folder=../trained_models/bmisner' %(tmp_path,tmp_path)
	os.system(outside_order)
	print(glob.glob(os.path.join(tmp_path,"tmp*")))
	result_path = glob.glob(os.path.join(tmp_path,"tmp*"))[0]
	result_file = os.path.join(result_path,'000_deploy.txt')
	re_tag = []
	with open(result_file,'r') as handle:
		for line in handle:
			a_line = line.strip()
			if a_line == '': continue
			a_line = a_line.split(" ")
			if(a_line[8]=="B-BME"): # new
				re_tag.append(a_line[0]) #new
#			re_tag.append({"term":a_line[0],"F1":a_line[4],"F2":a_line[5], "F3":a_line[6], "F4": a_line[7], "Result": a_line[8] })
		shutil.rmtree(tmp_path, ignore_errors=True)
		# return(json.dumps(re_tag)) #new 
		return(re_tag) #new 


# if __name__ == "__main__":
# 	kk=run_NER(sys.argv[1])
# 	print(kk)
# 	print(sys.argv[2])
# 	with open(sys.argv[2],"w") as fileout:
# 		for i in kk:
# 			fileout.write(i["term"]+"\t"+i["type"]+"\n")
