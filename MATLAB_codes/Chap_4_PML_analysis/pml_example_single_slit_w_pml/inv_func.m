function result=inv_func(R1,R2);

Iden=eye(size(R1,1));
result=inv(Iden-R1*R2);