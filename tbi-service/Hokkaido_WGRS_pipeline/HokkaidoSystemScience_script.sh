python HokkaidoSystemScience_script.py -f Raw/TN1802D0175/TN1802D0175_TKD180301122_HJCY7CCXY_L6_1.fq.gz -r Raw/TN1802D0175/TN1802D0175_TKD180301122_HJCY7CCXY_L6_2.fq.gz -s Ig13979
python HokkaidoSystemScience_script.py -f Raw/TN1802D0176/TN1802D0176_TKD180301123_HJCY7CCXY_L7_1.fq.gz -r Raw/TN1802D0176/TN1802D0176_TKD180301123_HJCY7CCXY_L7_2.fq.gz -s Ig13980
python HokkaidoSystemScience_script.py -f Raw/TN1802D0177/TN1802D0177_TKD180301124_HJ7KMCCXY_L3_1.fq.gz -r Raw/TN1802D0177/TN1802D0177_TKD180301124_HJ7KMCCXY_L3_2.fq.gz -s Ig13981
python HokkaidoSystemScience_script.2.py --finalBams Ig13979.final.bam Ig13980.final.bam Ig13981.final.bam
python HokkaidoSystemScience_script.3.py --ids Ig13979 Ig13980 Ig13981
