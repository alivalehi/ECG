clear all
load('training_data_tape222.mat');

data1 = linear_db.N;
data2 = linear_db.V;
data3 = linear_db.S;
data4 = linear_db.F;
[ ro1,ro2,ro3,ro4,Raw_sph,ortho_sph_cooked,CloseIndex,FarIndex] = train( data1,data2,data3,data4);
X = x_linear 
[ Y ] = predict( X,Raw_sph,ortho_sph_cooked,CloseIndex,FarIndex )
