function ret=find_best_fit(F,W,terms)
terms2=terms;
for j=1:size(terms,2)
	terms2(:,j)=terms2(:,j).*W;
end;
coeffs=inv(terms2'*terms2)*terms2'*(F.*W);
ret=terms*coeffs;