import time, datetime
from Bio import Entrez


def oneparser(pmid):
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



def localparser(filename):
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