

# coding: utf-8

# In[1]:


import pickle
from glob import glob
import os
import pandas as pd
import re
import argparse


def locate_rank(number, my_list):
    # Step 1: Sort the list in ascending order
    sorted_list = sorted(my_list)

    # Step 2: Insert the number into the sorted list while maintaining the sorted order
    from bisect import insort_left
    insort_left(sorted_list, number)

    # Step 3: Find the index of the inserted number in the list
    index = sorted_list.index(number)

    # Step 4: Calculate the percentage rank using the index and the length of the list
    percentage_rank = (index / len(sorted_list))

    return percentage_rank

def map_rank(input,score_dict):
    output=re.sub('no_rank','with_rank',input)
    print('the scores in: '+input+' is ranking to '+output)
    res_check = pd.read_csv(input)
    res_check['Rank'] = res_check.apply(lambda row: locate_rank(row['Score'], score_dict[row['Antigen']]), axis=1)
    res_check.to_csv(output,index=False)


if __name__ == "__main__":
    current = os.path.dirname(os.path.realpath(__file__))
    parent = os.path.dirname(current)


    parser = argparse.ArgumentParser(description='Parameters for adding ranks.')

    parser.add_argument('--background_score',type = str, help = 'the pkl file of the score dictionary of background BCRs',default = os.path.join(parent, 'data/background/background_score_dict_10000.pkl'))
    parser.add_argument('--input_to_rank',type = str,help = 'the results in csv with scores but no rank',default=os.path.join(parent, 'data/example/output/binding_results_no_rank.csv'))

    args = parser.parse_args()


    with open(args.background_score,'rb') as score:
        score_dict=pickle.load(score)

    ToRANK = args.input_to_rank
    map_rank(ToRANK,score_dict)

    #     res_check['Rank'].head()




# os.chdir('/project/DPDS/Wang_lab/shared/BCR_antigen/data/melanoma/')
#
# binding_files=glob('output_*/merged_results_no_rank.csv')
# print(len(binding_files))


# In[16]:


#binding_files


# In[17]:



# output='rank_results/'+'_'.join(re.sub('merged_results_no_rank.csv','binding_results_update.csv',binding_files[0]).split('/')[2:])


# In[18]:


# output


# In[19]:


#res_check.to_csv(output)


# In[20]:


# binding_files=glob('patch*/*/binding*csv')


# In[21]:


# print(len(glob('patch*csv')))
# print(len(binding_files))


# In[ ]:



# i = 0
# for binding in binding_files:
#     output=re.sub('merged_results_no_rank.csv','binding_results_update.csv',binding)
#     # print(output)
#     if os.path.exists(output):
#         print('the '+str(i)+'-th sample: '+binding+' has been ranked to '+output)
#     else:
#         print('the '+str(i)+'-th sample: '+binding+' is ranking to '+output)
#         res_check = pd.read_csv(binding,index_col=0)
#         res_check['Rank'] = res_check.apply(lambda row: locate_rank(row['Score'], score_dict[row['Antigen']]), axis=1)
# #     res_check['Rank'].head()
#         res_check.to_csv(output)
#     i += 1
