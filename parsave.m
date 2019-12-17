function parsave(PATH,saving_file_name, SNRdB,hChan,hEnc,hDec,BER,PER,delta)
save([PATH,saving_file_name,'.mat'],'SNRdB','hChan','hEnc','hDec','BER','PER','delta')
end

