import pandas as pd
import sys
from Bio import Entrez
Entrez.email = "seika.ooiwa@gmail.com"

#path情報
base_path =  '/Users/seika/Desktop/virtual_screening'

def create_query():
    """pubmed検索用のクエリー作成に用いるデータセット（機能名：定義）を生成
    
    Returns:
    --------
    category: list
        機能名のリスト
        
    func_name_to_query: dict
        機能名と定義（辞書型）ex. {"anti_bacteria/fungi":['antibacterial','antimicrobial','antifungal','antimycotic','fungicide'],,}
    """
    
    # 機能の定義
    func_name_to_keywords = {
        'anti_bacteria/fungi':['antibacterial','antimicrobial','antifungal','antimycotic','fungicide'],
        'insecticide':['insecticidal','insecticide','pesticide'],
        'herbicide':['herbicide','weed killer','weedkiller'],
        'tickkiller':['acaricidal','tick killer','tickkiller'],
        'virus':['virucidal','antivirus','antiviral'],
    }
    # queryの生成
    func_name_to_query = {
        func_name: get_query(keywords) for func_name, keywords in func_name_to_keywords.items()
    }

    category = list(func_name_to_query.keys())

    return category,func_name_to_query
    
def get_query(keywords,field="[Title/Abstract]"):
    """keywords中の機能名リストを取り出し、field名と連結した上で、機能名同しをOR連結する
    
    Parameters
    ----------
    keywords: list
        機能名のリスト
    field: str
        Pubmed上の検索領域（default: Title/Abstract)
    Returns
    -------
    "(" + ors + ")": str
        検索文字列
    """
    ors = " OR ".join(f"({keyword}{field})" for keyword in keywords)
    return "(" + ors + ")"

def check_spell(tname):
    """
    """
    kyw = Entrez.read(Entrez.espell(term=tname))
    kyw2 = kyw['CorrectedQuery']

    return kyw2

def search_articles(search_word, max_results=5):
    
    """keywordを元にPubmedを検索しヒットしたpubidを取得

    Parameters
    ----------
    serch_word: str
        検索キーワード("target_nane AND [func1 OR func2]")       
    max_results: int
        検索結果の表示数, デフォルト(5件)

    Returns
    -------
    results: dict
    """
    
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax=str(max_results),
                            term=search_word)
    results = Entrez.read(handle,validate=False)
    return results

def fetch_abstract(result):

    """pubidを元にAbstractを読み込み

    Parameters
    ----------
    result: dict
        EntrezでPubmedのキーワード検索結果 
    
    Returns
    -------
    dataframe: dataframe
    """
    
    data = []
    
    for i in range(len(result['IdList'])):
        pmid = result['IdList'][i]
    
        handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
        article = Entrez.read(handle)['PubmedArticle'][0]

        try:
            abstract = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            add = [pmid,abstract]
            data.append(add)

    
        except KeyError as e:
            print('keyError')
            
    dataframe = pd.DataFrame(data,columns=['pubid','Abstract'])
    
    return dataframe

def calculate_data(abdf,keys):

    #abdfが空データフレームか判定
    if abdf.empty:
        return '-',abdf

    else:
        #keysを個々の検索ワードに分解
        #keys2 = keys.replace('AND',',').replace('OR',',').replace(' ','').split(',')
        keys2 = [word.strip('-.()/[]') for word in keys.replace('AND',',').replace('OR',',').replace(' ','').replace('[Title/Abstract]','').split(',')]
        #dataframeにカウント結果入力カラム追加
        abdf[keys2]=''
            
        #abstract中のキーワードの出現回数をカウントし、dataframeに入力
        for i in range(len(abdf)):
            abst_text = [word.strip('-.()/') for word in abdf.loc[i,'Abstract'].lower().split()]
            for k in range(len(keys2)):
                sword = keys2[k]
                ct = abst_text.count(sword)
                abdf.loc[i,sword] = int(ct)
    
        
        #検索ワードを候補物質名と機能リストに分解
        target_1 = keys2.pop(0)
        target_2 = keys2.pop(0)
        func_list2 = keys2
        
        #候補物質名のみや機能名のみでヒットしたabstractを削除
        tmp = abdf[func_list2]
        abdf['sum_func'] = tmp.sum(axis=1)
        abdf2 = abdf[(abdf['sum_func']!=0)]
        abdf3 = abdf2[(abdf2[target_1]!=0) | (abdf2[target_2]!=0)].reset_index()
        
               
        #候補物質が機能名に関係するか確認
        #正確性を上げるため、２件以上の論文があることを条件にする
        if len(abdf3) >= 2:
            return 'Hit',abdf
    
        else:
            return '-',abdf

def search_metabolites_in_pubmed(asdf):
    """dataframeに記載の代謝物名と機能をクエリ―としてpubmed検索を行い、dataframeに結果を追加する
    
    Parameter:
    ----------
    asdf: dataframe
        列名='Candidate_metabolite'の下に代謝物名が記載
    
    Return:
    -------
    asdf: dataframe
    raw_data: dataframe
    
    See Also:
    --------
    get_query(keywords,field="[Title/Abstract]")
    check_spell(tname)
    search_articles(serch_word, max_results=5)
    fetch_abstract(result)
    calculate_data(abdf,keys)
    
    """
    # クエリ（機能）の生成
    category,func_name_to_query = create_query()
    
    # データフレームに機能名のカラムを追加
    asdf[category]=''
    
    raw_data = []

    for i in range(len(asdf)):
        # 検索のため小文字に統一
        tname = asdf.iloc[i]['Candidate_metabolite'].lower()
        # spellcheck
        tname_s = check_spell(tname)

    #　関連文献のabstract抽出
        for func_category in category:

            #検索式の生成
            func_word = func_name_to_query[func_category]
            search_word = f'(({tname_s}[Title/Abstract]) OR ({tname}[Title/Abstract])) AND {func_word}'
            #pubmed検索
            result = search_articles(search_word)
            #abstractの抽出(dataframe)
            abs_df = fetch_abstract(result)
            #候補物質と機能の関係性評価(rs;<1:関連あり,0:関連なし>)
            rs,abs_df = calculate_data(abs_df,search_word)
            #Antismash解析データフレームに入力
            asdf.loc[i,func_category]=rs
            #row_dataを保存
            raw_data.append(abs_df)

    return asdf,raw_data 

file_path = sys.argv[1]
output_path = sys.argv[2]

df = pd.read_csv(f'{file_path}')
df,raw_data = search_metabolites_in_pubmed(df)

df.to_csv(output_path)

