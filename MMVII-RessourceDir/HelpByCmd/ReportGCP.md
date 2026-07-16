# Help for command ReportGCP

**Table of Contents**
- [Description](#description)
- [VisuResult](#visualresult)
- [StatisticalResult](#statisticalresult)


# Description
This command generate report on quality of orientation relative to a set of Ground 
control point. It generate two kind of information :
- ** statisical information in csv format , these file are generated in the  ReportDir of the application
- ** visual information that contain a visualization of the residual for each camera an, optionnaly, for each image,
     these result are stored in the VisuDir of the application

#statisticalresult

Generate the following csv file
- ** ByCam.csv  : residual of projection in pixel, agregated by camera, usable to analyse internal calibration quality
- ** ByGCP.csv  : residual of projection in pixel, agregated by 3D point, usable to analyse internal calibration quality
- ** ByGCP3D.csv  : differnce in ground coordinate between 3D coordinate and budle intesection
- ** ByGCP3D\_Stat.csv  : statistique on valude of  ByGCP3D.csv (average ...)
- ** ByImage.csv : residual of projection in pixel, agregated by images, usable to analyse pose estimation quality
- ** Detail.csv : residual of projection in pixel for each point/image, usable to analyse outliers


#visualresult
For each camera ABC, the following low resolution images are generated to anlyse the distribution of GCP on the sensor and 
the bias of reprojection :

- ** W\_RawResidual\_ABC  : raw distribution of point on the sensor plane (generally sparse )
- ** W\_FiltResidual\_ABC  : filtered distribution of point on the sensor plane (smooth)
- **  X\_Residual\_ABC and  Y\_Residual\_ABC average of reprojection error, usable to detect spatial bias

When the optional vector GFV (generate vector field) is used, a  visualaition  is genrated for each image ABC, that show for each GCP if it is
detected and, if it is,  show the vector between detection and projection . 

- ** note GFV= [Mul,Witdh,Ray,Zoom?=2,JPeg?=1]
- ** for each detected GCP, a vector is drawn with an amplification of Mul and a witdh Witdh
- ** for each undected GGP, a circle is drawn a the center of projection  and  of ray Ray
- ** Zoom indicate the reduction factor
- ** JPeg indicate if Jpeg/Tiff format is used





