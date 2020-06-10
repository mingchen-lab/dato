import pymysql
import struct
import socket
import pickle
import re
import io
import urltools
import time, datetime

def url_extractor(text):
    http_or_https = re.compile(r"(?isu)(https?\://[a-z0-9\.\?/&\=\:\#\-_]+)")
    ftp = re.compile(r"(?isu)(ftp\://[a-z0-9\.\?/&\=\:\#\-_]+)")
    wrong_http = re.compile(r"(?isu)(http\:?@[a-z0-9\.\?/&\=\:\#\-_]+)")
    wrong_ftp = re.compile(r"(?isu)(ftp\:@[a-z0-9\.\?/&\=\:\#\-_]+)")
    http_url = http_or_https.findall(text)
    ftp_url = ftp.findall(text)
    http_with_special = wrong_http.findall(text)
    ftp_with_special = wrong_ftp.findall(text)
    http_recorrect = []
    for a_wrong_http in http_with_special:
        real_wrong_http = a_wrong_http.rstrip('.')
        if real_wrong_http[4] == '@':
            real_http = r"http://" + real_wrong_http[5:]
            http_recorrect.append(real_http)
        if real_wrong_http[4:6] == ":@":
            real_http = r"http://" + real_wrong_http[6:]
            http_recorrect.append(real_http)
        else:
            pass
    ftp_recorrect = []
    for a_wrong_ftp in ftp_with_special:
        real_wrong_ftp = a_wrong_ftp.rstrip('.')
        if real_wrong_ftp[3] == '@':
            real_ftp = r"ftp://" + real_wrong_ftp[5:]
            ftp_recorrect.append(real_ftp)
        if real_wrong_ftp[3:5] == ":@":
            real_ftp = r"ftp://" + real_wrong_ftp[6:]
            ftp_recorrect.append(real_ftp)
        else:
            pass
    result = (http_url+ftp_url+http_recorrect+ftp_recorrect)
    return([urltools.normalize(r_norm).rstrip("\.") for r_norm in result])

def Urltools_dict(url):
    arser_urltools = urltools.parse(url)
    parser_dict = {}
    parser_dict["scheme"]=arser_urltools.scheme
    parser_dict["username"]=arser_urltools.username
    parser_dict["password"]=arser_urltools.password
    parser_dict["subdomain"]=arser_urltools.subdomain
    parser_dict["domain"]=arser_urltools.domain
    parser_dict["tld"]=arser_urltools.tld
    parser_dict["port"]=arser_urltools.port
    parser_dict["path"]=arser_urltools.path
    parser_dict["query"]=arser_urltools.query
    parser_dict["fragment"]=arser_urltools.fragment
    parser_dict["url"]=arser_urltools.url
    return(parser_dict)


def domain2ip10(hosts):
    ip = socket.gethostbyname(hosts)
    packedIP = socket.inet_aton(ip)
    return(struct.unpack("!L", packedIP)[0])

import pymysql
conn = pymysql.connect(host="localhost",user="root",passwd="*",db="ip2geo")
cur = conn.cursor()

def ip2location(ip):
    conn = pymysql.connect(host="localhost",user="root",passwd="fox12345",db="ip2geo")
    cur = conn.cursor()
    sql = "select * from ip2location_db11 where ip_to >= %s and ip_from <= %s limit 1"
    cur.execute(sql,(ip,ip))
    return(cur.fetchone())

def meta_url(hosts):
    r = {}
    r['url'] = hosts
    url_dict = Urltools_dict(hosts)
    r['url_meta'] = url_dict
    if len(hosts.split("zju.edu"))>1:
        r["urlinfo"]={"country":"China","city":"Hangzhou","state":"Zhejiang","urlloc":{"lat":"30.29365","lng":"120.16142"}}
        r['country']="CN"
        r['urlloc']={"lat":"30.29365","lng":"120.16142"}
        return(r)
    try:
        if url_dict['subdomain']!='':
            domain = url_dict['domain'] + "." +url_dict['tld']
            subdomain = url_dict['subdomain'] + "." +  url_dict['domain'] + "." + url_dict['tld']
        else:
            domain = url_dict['domain'] + "." +url_dict['tld']
            subdomain = url_dict['domain'] + "." +url_dict['tld']
        r['or_domain_info'] = ip2location(domain2ip10(subdomain))
        r['domain_info'] = ip2location(domain2ip10(subdomain))
        r["urlinfo"]={"country":r['domain_info'][3],"city":r['domain_info'][5],"state":r['domain_info'][4],"urlloc":{"lat":r['domain_info'][6],"lng":r['domain_info'][7]}}
        r['country']=r['domain_info'][3]
        r['urlloc']={"lat":r['domain_info'][6],"lng":r['domain_info'][7]}
    except:
        r['or_domain_info'] = None
        r['domain_info'] = None
    return(r)
