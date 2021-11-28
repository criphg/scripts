import urllib.request
import urllib.parse

idx_v = 141255

def get_sequence_addgene(id_v):
    url = 'http://www.addgene.org/'+str(id_v)+'/sequences/'

    headers = {}
    headers['User-Agent'] = "Mozilla/5.0"

    req = urllib.request.Request(url, headers=headers)
    resp = urllib.request.urlopen(req)
    #resp_data = resp.read()

    seq = ""
    tag_v=0
    for line in resp:
        if str(line).find(">&gt; Addgene NGS Result") == -1:
            continue
        else:
            while(str(line).find("</textarea>") == -1):
                if(line != ""):
                    decoded_line = line.decode("utf-8")
                    #decoded_line = line
                    #print(decoded_line)
                    line = resp.readline()
                    if(tag_v == 0):
                        seq = ">"+str(id_v)+"\n"

                    else:
                        seq = seq + decoded_line.strip()
                    tag_v = 1

    #print(resp_data)
    return(seq)
    

fasta_seq = get_sequence_addgene(idx_v)
print(fasta_seq)