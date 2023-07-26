import pandas

metadata = pandas.read_csv('data/metadata.csv',
                           sep=',',
                           na_values=''
                           )

metadata.rename(columns={'Collection_Date':'date'}, inplace=True)

metadata.replace(float('nan'),'?',inplace=True)

division_column=[]
for i in metadata.Geo_Location:
    #geo location will become region
    if i != '?':
        i=i.split(':')
        if len(i)>1:
            division=i[1]
            if division[0]==' ':
                division=division[1:]
            e=division
        else:
            e='?'
    division_column.append(e)

metadata['Division']=division_column
metadata.drop('Geo_Location',axis=1)

metadata.drop_duplicates(subset='Accession')

metadata.to_csv('data/processed_metadata.csv')