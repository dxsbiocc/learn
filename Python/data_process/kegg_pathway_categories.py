import warnings
import requests
from lxml.etree import HTML


def get_pathways(url_pathways, url_links):
    # get pathway list page
    req = requests.get(url_pathways)
    htree = HTML(req.text)
    # extract pathway info
    hsa = []
    for e in htree.cssselect('#{p}, #{p} ~ b, #{p} ~ h4'.format(p="metabolism")):
        if e.tag == 'h4':
            cate = e.text.split(' ', maxsplit=1)[1]
        elif e.tag == 'b':
            subcate = e.text.split(' ', maxsplit=1)[1].split('-')[0]
            div = e.getnext()
            for t in div.cssselect('dt, a'):
                if t.tag == 'dt':
                    hsa.append([cate.strip(), subcate.strip(), f'hsa{t.text.strip()}', None])
                elif t.tag == 'a':
                    if hsa[-1][-1]: warnings.warn('Pathway has been exists!', UserWarning)
                    hsa[-1][-1] = t.text.strip()
        else:
            raise ValueError('Unexpected value of', e.tag)
    # to DataFrame
    data = pd.DataFrame(hsa, columns=['categories', 'subcategories', 'ID', 'name'])
    # get links between pathway and gene
    hsa = pd.read_csv(url_links, sep='\t', names=['gene', 'path'])
    hsa = hsa.applymap(lambda x: x.split(':')[1])
    # remove prefix and rearange the links
    tmp = hsa.groupby('path').apply(lambda x: pd.Series([x.name, ','.join(x['gene']), len(x['gene'])], 
                                                        index=['ID', 'gene_list', 'count'])).reset_index(drop=True)
    # merge
    pathway = pd.merge(data, tmp)
    return pathway

if __name__ == '__main__':
    url_pathways = 'https://www.kegg.jp/kegg/pathway.html'
    url_links = 'https://rest.kegg.jp/link/pathway/hsa'
    pathway = get_pathways(url_pathways, url_links)
