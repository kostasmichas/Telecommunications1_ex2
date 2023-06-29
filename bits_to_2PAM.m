function b2p = bits_to_2PAM(b)
  b2p=zeros(length(b),1);
  for i=1:length(b)
    if(b(i)==0)
      b2p(i)=1;
    elseif(b(i)==1)
      b2p(i)=-1;
    end
  end
end