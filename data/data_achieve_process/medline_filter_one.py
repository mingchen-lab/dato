import time, datetime
from Bio import Entrez
from medline import parse
import re
import pickle

def Url_nonsense(url):
	nonsense_words = ["ow.ly","lww.com","clinicaltrial","irct.ir","aacrjournals.org","doi.org","creativecommons.org",
					  "trialregister.nl","webcitation.org","aacrjournals.org","surgicalneurologyint.com","annualreviews.org",
					  "wiley.com","elseviercme.com","biomedcentral.com","sagepub.com","rimed.org","elsevier.com","link.springer-ny.com",
                      "www.biologists.com","anzctr.org.au","www.aseronline.org","w3.org","chictr.org","bit.ly","chinadrugtrials","ncdirindia",
                      "umin.ac.jp","isrctn.com","wileyhealthlearning","disease-perception","fallsclinic"]
	for word in nonsense_words:
		try:
			if re.search(word,url,re.I).span():
				return(True)
			else:
				continue
		except:
			pass
	return(False)


def TI_nonsense(TI):
    nonsense_words = ["WITHDRAWN","RETRACTED:","Re:","Core curriculum illustration","Erratum:","Retracted:","Retraction",
                      "Educational Case:","REMOVED", "TEMPORARY REMOVAL","ftp","Rejected","\[","\?","Commentary","case report",
                     "clinical trial","clinical trial","case study","survey"," review ","official publication","comparison study"
                     "Annual Report","phase I","official journal","Clinical & experimental metastasis","Health Clinical","Cerebral cortex",
                     "congenital heart disease","preliminary study","a review","controlled trial","Clinical use","summary of",
                     "technical note on","Structural analysis","who"," what ","surgery","a randomised crossover trial","etrospective report","a systematic review"]
    for word in nonsense_words:
        try:
            if bool(re.search(word,TI,re.I).span()):
                return(True)
            else:
                pass
        except:
            pass
    return(False)

def JT_nonsense(JT):
    nonsense_words = ["Plant disease","Experimental and clinical immunogenetics","Radiology","Nature medicine","Cerebral cortex",
                     "The Journal of clinical endocrinology and metabolism","Diagnostic pathology","The Journal of clinical investigation",
                     "Global change biology","Trials","Neurosurgical focus","Current opinion in investigational","Fibrogenesis & tissue repair",
                     "Chest","European journal of nutrition","Indian journal of dermatology","Pilot and feasibility studies","Molecular plant pathology",
                     "Circulation research","Journal of the American Heart Association","Journal of burn care & research","Blood pressure monitoring",
                     "Circulation","Genes & cancer","Water science and technology","Advances in experimental medicine and biology","Acta naturae",
                     "Acta histochemica","Acta crystallographica. Section D, Structural biology","ACS synthetic biology","ACS nano","Acoustics Australia",
                     "Accounts of chemical research","3 Biotech","AIDS and behavior","American journal of physical medicine & rehabilitation","American journal of physiology. Heart and circulatory physiology",
                     "Aesthetic surgery journal","Age and ageing","Age","Aging","AIDS and behavior","AIDS and behavior","AIDS research and human retroviruses",
                     "AJR","Alimentary pharmacology & therapeutics","Allergy","Alzheimer","AMIA Symposium","Anaesthesia","Anais da Academia","Analytical cellular pathology",
                     "Anatomical record","Atherosclerosis","Biological","Biology of","Biology open","Blood","BMC anesthesiology","BMC musculoskeletal disorders","BMC neurology",
                     "BMJ open","BMJ quality improvement","BMJ supportive & palliative care","Bone","Breast","British","Canadian",
                      "Cancer biology & therapy","Chinese medical journal","Clinica","Critical care","Diabetic medicine","Diabetologia","Diagnosis",
                      "Epilepsia"," surgery","FASEB journal","The Journal of the Acoustical Society of America","The Journal of urology","The oncologist",
                     "The Review of scientific instruments","The Science of the total environment"]
    for word in nonsense_words:
        try:
            if bool(re.search(word,JT,re.I).span()):
                return(True)
            else:
                pass
        except:
            pass
    return(False)

def oneparser(pmid):
    addic = pickle.load(open("address_version3_7.pkl","rb"))
    handle=Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
    res = parse(handle)
    out = []
    for r in res:
        if "AD" in r.keys():
            r["MAIL"]=EMmail_r(r['AD'])
            r["AD"] = raw_medline_ad_washer(r['AD'])
        if 'CRDT' in r.keys():
            conv = time.strptime( r['CRDT'][0], "%Y/%m/%d %H:%M" )
            r['CRDT'] = datetime.datetime(*conv[:6]) # date created
        if "PMID" in r.keys()  and r["PMID"] != "":
            r['_id'] = int(r['PMID'])
            try:
                out.append(r)
            except Exception as e:
                print(e)
                continue
        else: 
            continue
    return(out)

def multi_parser(pmidlist):
    addic = pickle.load(open("address_version3_7.pkl","rb"))
    handle=Entrez.efetch(db="pubmed", id=pmidlist, rettype="medline", retmode="text")
    res = parse(handle)
    out = []
    for r in res:
        if "AD" in r.keys():
            r["MAIL"]=EMmail_r(r['AD'])
            r["AD"] = raw_medline_ad_washer(r['AD'])
        if 'CRDT' in r.keys():
            conv = time.strptime( r['CRDT'][0], "%Y/%m/%d %H:%M" )
            r['CRDT'] = datetime.datetime(*conv[:6]) # date created
        if "PMID" in r.keys()  and r["PMID"] != "":
            r['_id'] = int(r['PMID'])
            try:
                out.append(r)
            except Exception as e:
                print(e)
                continue
        else: 
            continue
    return(out)

def localparser(filename):
    addic = pickle.load(open("address_version3_7.pkl","rb"))
    handle=open(filename,'r')
    res = parse(handle)
    out = []
    for r in res:
        if "AD" in r.keys():
            r["MAIL"]=EMmail_r(r['AD'])
            r["AD"] = raw_medline_ad_washer(r['AD'])
        if 'CRDT' in r.keys():
            conv = time.strptime( r['CRDT'][0], "%Y/%m/%d %H:%M" )
            r['CRDT'] = datetime.datetime(*conv[:6]) # date created
        if "PMID" in r.keys()  and r["PMID"] != "":
            r['_id'] = int(r['PMID'])
            try:
                out.append(r)
            except Exception as e:
                print(e)
                continue
        else: 
            continue
    return(out)



def raw_medline_ad_washer(AD):
    addic = pickle.load(open("address_version3_7.pkl","rb"))
    # remove email
    mail_word = [". Electronic address:","Contact:contact-",
                 "paragraph sign paragraph sign",". CONTACT:",". Email:",". email:",". E-mail:","; Electronic address:","; CONTACT:","; Email:","; email:","; E-mail:"]
    for word in mail_word:
        AD=AD.replace(word," ")
    AD = AD.replace("  "," ")
    AD = re.sub(r'((\d){4}\s(\d){4}\s(\d){4}\s(\d){4})','',AD)
    AD = re.sub(r'grid\.\S+','',AD)
    AD = re.sub(r'(\d){14}\S+','',AD)
    AD = re.sub(r'(\w\.){3,}','',AD)
    AD = re.sub(r'\[(\d){1,2}\]','',AD)
    AD = re.sub(r'\[(\d){1,2}\]','',AD)
    AD = re.sub(r'(\d){1,2}\]','',AD)
    email = EMmail_r(AD)
    for i in email:
        AD = re.sub(i,"",AD)
    # remove (***)
    AD = re.sub('(\s){1,}'," ",AD)
    AD = re.sub('\((\s|\w|\d|\*|\.)+\)',"",AD)
    for i in addic['coun_non']:
        AD = AD.replace(i,i.replace('.',"***")) 
    for i in addic['coun_and']:
        AD = AD.replace(i,i.replace(' and',"."))

    AD = AD.replace("  "," ").replace("; .",".").replace(". .",'.').replace("..",".").replace("<>","").replace(". ;",'.').replace(", ,",",").replace(",,",",").replace("()","")
    AD = AD.replace(';','. @@@').replace('. @@@',".##").replace(' @@@','.##')   
    res = [k.strip().rstrip(".").split(".## " or ';' )  for k in AD.strip().rstrip(".").split(". " or "")]
    res1=[]
    for i in res:
        res1.append([re.sub('\S+@\S+',"",n,flags=re.M).replace("***",".").replace(" ,",",").replace("###### ","").strip().rstrip(".") for n in i])
    return(res1)

def EMmail_r(AD):
    return(list(set(i.strip().rstrip(".") for i in re.findall('\S+@\S+',AD.replace("@@@","###").replace("***",".")))))