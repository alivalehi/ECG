function output = PCA_visualization( test_data,D )

[wcoeff,score,latent,tsquared,explained] = pca(test_data(:,1:D));
data_3D = test_data(:,1:D)*wcoeff(:,1:3);

figure
scatter3(data_3D(:,1),data_3D(:,2),data_3D(:,3),[],test_data(:,D+1))
output = [data_3D(:,1) data_3D(:,2) data_3D(:,3) test_data(:,D+1)];
end

