clear all
load('training_data_tape222.mat');

data1 = linear_db.N;
data2 = linear_db.V;
data3 = linear_db.S;
data4 = linear_db.F;
nd1= size(data1,1);
nd2= size(data2,1);
nd3= size(data3,1);
nd4= size(data4,1);
Fisher( [data1',data2',data3',data4'],[1,nd1,nd2+nd1,nd3+nd2+nd1,nd4+nd3+nd2+nd1])
[ ro1,ro2,ro3,ro4,Raw_sph,ortho_sph_cooked,CloseIndex,FarIndex] = train( data1,data2,data3,data4);
nd1= size(ro1,1);
nd2= size(ro2,1);
nd3= size(ro3,1);
nd4= size(ro4,1);
Fisher( [ro1',ro2',ro3',ro4'],[1,nd1,nd2+nd1,nd3+nd2+nd1,nd4+nd3+nd2+nd1])
    acosd(1-pdist2(nsphere2cart(ortho_sph_cooked(2,:)),nsphere2cart(ortho_sph_cooked(1,:)),'cosine'));
    acosd(1-pdist2(nsphere2cart(ortho_sph_cooked(3,:)),nsphere2cart(ortho_sph_cooked(1,:)),'cosine'));
X = x_linear; 
[ Y ] = predict( X,Raw_sph,ortho_sph_cooked,CloseIndex,FarIndex )
