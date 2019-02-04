load('WL_sampl.dat');
load('WL_exh.dat');
load('WL_5x5.dat');
load('WL_grid.dat');
load('Grid_UC.dat');

x=WL_sampl;
G=WL_grid;
R=WL_5x5;
G_UC=Grid_UC;
block=[5 5];
nd=[5 5];
ival=0;
nk=3;
nk_uc=120;
rad=99999;
ntok=1;
avg=mean(x(:,3));
smu=[5 5];
panneau=[65 75];
vc=[0 0.05:9:1050];

for cas=3:3
    cas
    switch cas
        case 1
            model=[1 1;2 15];
            c=[2000;80000];
            
        case 2
            model=[1 1;2 10];
            c=[16000;74000];
            
        case 3
            model=[1 1;2 12];
            c=[5000;85000];
    end
    
    [stat,x0s_ok,x0s_ck,ton_uc,mean_uc]=okckuc(x,G,R,model,c,block,nd,ival,nk,nk_uc,rad,ntok,avg,smu,panneau,vc,G_UC,cas);
    
    s_OK=latex(vpa(sym(stat.OK),5))
    s_CK=latex(vpa(sym(stat.CK),5))
    
end





