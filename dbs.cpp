bool pixkit::halftoning::iterative::dualmetricDBS2002(const cv::Mat &src1b, cv::Mat &dst1b){

	//////////////////////////////////////////////////////////////////////////
	/// exceptions
	if(src1b.type()!=CV_8UC1){
		assert(false);
	}

	//////////////////////////////////////////////////////////////////////////
	///// get filter
	vector<Mat>	cpp1d;
	cpp1d.push_back(Mat());
	cpp1d.push_back(Mat());
	pixkit::halftoning::ungrouped::generateTwoComponentGaussianModel(cpp1d[0],43.2,38.7,0.0219,0.0598);
	pixkit::halftoning::ungrouped::generateTwoComponentGaussianModel(cpp1d[1],19.1,42.7,0.0330,0.0569);
	int FilterSize	=	cpp1d[0].rows;
	int	exFS=FilterSize;
	int	tempFS=FilterSize/2;

	//////////////////////////////////////////////////////////////////////////
	///// get weight
	Mat	weightmap1f(Size(2,256),CV_32FC1);
	for(int i=0;i<256;i++){
		float	gray	=	(float)i/255.;
		// get weight
		if(gray<=0.25){
			weightmap1f.ptr<float>(i)[0]	=	std::sqrtf(1-((float)4.*gray-1.)*((float)4.*gray-1.));
		}else if(gray>0.25&&gray<=0.75){
			weightmap1f.ptr<float>(i)[0]	=	std::fabsf((float)4.*gray-2);
		}else{
			weightmap1f.ptr<float>(i)[0]	=	std::sqrtf(1-((float)4.*gray-3.)*((float)4.*gray-3.));
		}
		weightmap1f.ptr<float>(i)[1]	=	1.-weightmap1f.ptr<float>(i)[0];
	}

	//////////////////////////////////////////////////////////////////////////
	/// initialization
	int	m_Height	=	src1b.rows;
	int	m_Width		=	src1b.cols;
	Mat	dst1f;
	dst1f.create(src1b.size(),CV_32FC1);

	//////////////////////////////////////////////////////////////////////////
	// get halftone image
	srand(0);
	for(int i=0;i<m_Height;i++){
		for(int j=0;j<m_Width;j++){
			double	temp=((double)rand())/32767.;
			dst1f.ptr<float>(i)[j]=temp<0.5?0.:255.;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	/// Change grayscale to absorb
	for(int i=0;i<m_Height;i++){
		for(int j=0;j<m_Width;j++){
			dst1f.ptr<float>(i)[j]	=	1.-dst1f.ptr<float>(i)[j]/255.;
		}
	}
	/// get error matrix
	Mat	em1f(Size(m_Width,m_Height),CV_32FC1);
	for(int i=0;i<m_Height;i++){
		for(int j=0;j<m_Width;j++){
			double	oriv			=	1.-((double)src1b.ptr<uchar>(i)[j])/255.;
			em1f.ptr<float>(i)[j]=dst1f.ptr<float>(i)[j]-oriv;
		}
	}
	/// get cross correlation
	vector<Mat>	c_ep1d;
	c_ep1d.push_back(Mat(Size(m_Width,m_Height),CV_64FC1));
	c_ep1d.push_back(Mat(Size(m_Width,m_Height),CV_64FC1));
	c_ep1d[0].setTo(0);
	c_ep1d[1].setTo(0);
	for(int i=0;i<m_Height;i++){
		for(int j=0;j<m_Width;j++){
			for(int m=i-tempFS;m<=i+tempFS;m++){
				for(int n=j-tempFS;n<=j+tempFS;n++){
					if(m>=0&&m<m_Height&&n>=0&&n<m_Width){
						c_ep1d[0].ptr<double>(i)[j]+=em1f.ptr<float>(m)[n]*cpp1d[0].ptr<double>(tempFS+m-i)[tempFS+n-j];
						c_ep1d[1].ptr<double>(i)[j]+=em1f.ptr<float>(m)[n]*cpp1d[1].ptr<double>(tempFS+m-i)[tempFS+n-j];
					}
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	///// DBS process
	int		BenefitPixelNumber;
	double	delta_E[10],a0[10],a1[10];
	while(1){
		BenefitPixelNumber=0;
		for(int i=0;i<m_Height;i++){	// entire image
			for(int j=0;j<m_Width;j++){

				//////////////////////////////////////////////////////////////////////////
				///// get weight
				int		currv	=	cvRound((float)dst1f.ptr<float>(i)[j]*255.);

				// = = = = = trial part = = = = = //
				// initialize psnr		0: original psnr, 1~8: Swap, 9: Toggle.
				// 8 1 2
				// 7 0 3
				// 6 5 4
				for(int m=0;m<10;m++){
					delta_E[m]=0.;
					a0[m]=0.;
					a1[m]=0.;
				}
				// change the delta error as per different replacement methods
				for(int mode=1;mode<10;mode++){
					int		m,n;
					if(mode>=1&&mode<=8){
						// set position
						if(mode==1){
							m=1;	n=0;
						}else if(mode==2){
							m=1;	n=1;
						}else if(mode==3){
							m=0;	n=1;
						}else if(mode==4){
							m=-1;	n=1;
						}else if(mode==5){
							m=-1;	n=0;
						}else if(mode==6){
							m=-1;	n=-1;
						}else if(mode==7){
							m=0;	n=-1;
						}else if(mode==8){
							m=1;	n=-1;
						}
						// get dE
						if(i+m>=0&&i+m<m_Height&&j+n>=0&&j+n<m_Width){
							// get weight to neighbor
							int	currv_nei	=	cvRound((float)dst1f.ptr<float>(i+m)[j+n]*255.);

							// get error
							if(dst1f.ptr<float>(i)[j]==1){
								a0[mode]=-1;
							}else{
								a0[mode]=1;
							}
							if(dst1f.ptr<float>(i+m)[j+n]==1){
								a1[mode]=-1;
							}else{
								a1[mode]=1;
							}
							if(dst1f.ptr<float>(i)[j]!=dst1f.ptr<float>(i+m)[j+n]){
								for(int w_idx=0;w_idx<2;w_idx++){
									delta_E[mode]+=(a0[mode]*a0[mode]	*	weightmap1f.ptr<float>(currv)[w_idx]*weightmap1f.ptr<float>(currv)[w_idx]	+	a1[mode]*a1[mode]	*	weightmap1f.ptr<float>(currv_nei)[w_idx]*weightmap1f.ptr<float>(currv_nei)[w_idx])	*	cpp1d[w_idx].ptr<double>(tempFS)[tempFS]	+
										2.*	a0[mode]	*	weightmap1f.ptr<float>(currv)[w_idx]	*	a1[mode]	*weightmap1f.ptr<float>(currv_nei)[w_idx]	*	cpp1d[w_idx].ptr<double>(tempFS+m)[tempFS+n]	+
										2.*	a0[mode]	*	weightmap1f.ptr<float>(currv)[w_idx]	*	c_ep1d[w_idx].ptr<double>(i)[j]	+
										2.*	a1[mode]	*	weightmap1f.ptr<float>(currv_nei)[w_idx]	*	c_ep1d[w_idx].ptr<double>(i+m)[j+n];
								}
							}
						}
					}else if(mode==9){
						if(dst1f.ptr<float>(i)[j]==1){
							a0[mode]=-1;
						}else{
							a0[mode]=1;
						}

						for(int w_idx=0;w_idx<2;w_idx++){
							delta_E[mode]	+=	(a0[mode]*a0[mode]	*	weightmap1f.ptr<float>(currv)[w_idx]*weightmap1f.ptr<float>(currv)[w_idx]	)	*	cpp1d[w_idx].ptr<double>(tempFS)[tempFS]	+
								2.*	a0[mode]	*	weightmap1f.ptr<float>(currv)[w_idx]	*	c_ep1d[w_idx].ptr<double>(i)[j];
						}

					}
				}
				// get minimum delta error and its position
				int		tempMinNumber	=0;
				double	tempMindE		=delta_E[0];
				for(int x=1;x<10;x++){
					if(delta_E[x]<tempMindE){
						tempMindE		=delta_E[x];
						tempMinNumber	=x;
					}
				}

				// = = = = = update part = = = = = //
				if(tempMindE<0){	// error is reduce
					// update hft image
					dst1f.ptr<float>(i)[j]	=1.-dst1f.ptr<float>(i)[j];
					if(tempMinNumber>=1&&tempMinNumber<=8){
						// get position
						int nm,nn;
						if(tempMinNumber==1){
							nm=1;	nn=0;
						}else if(tempMinNumber==2){
							nm=1;	nn=1;
						}else if(tempMinNumber==3){
							nm=0;	nn=1;
						}else if(tempMinNumber==4){
							nm=-1;	nn=1;
						}else if(tempMinNumber==5){
							nm=-1;	nn=0;
						}else if(tempMinNumber==6){
							nm=-1;	nn=-1;
						}else if(tempMinNumber==7){
							nm=0;	nn=-1;
						}else if(tempMinNumber==8){
							nm=1;	nn=-1;
						}
						// update hft image
						dst1f.ptr<float>(i+nm)[j+nn]	=1.-dst1f.ptr<float>(i+nm)[j+nn];
						// get weight to neighbor
						int	currv_nei	=	cvRound((float)dst1f.ptr<float>(i+nm)[j+nn]*255.);

						// update cross correlation
						for(int m=-tempFS;m<=tempFS;m++){
							for(int n=-tempFS;n<=tempFS;n++){
								if(i+m>=0&&i+m<m_Height&&j+n>=0&&j+n<m_Width){
									for(int w_idx=0;w_idx<2;w_idx++){
										c_ep1d[w_idx].ptr<double>(i+m)[j+n]+=a0[tempMinNumber]*weightmap1f.ptr<float>(currv)[w_idx]*cpp1d[w_idx].ptr<double>(tempFS+m)[tempFS+n];
									}
								}
								if(i+m+nm>=0&&i+m+nm<m_Height&&j+n+nn>=0&&j+n+nn<m_Width){
									for(int w_idx=0;w_idx<2;w_idx++){
										c_ep1d[w_idx].ptr<double>(i+m+nm)[j+n+nn]+=a1[tempMinNumber]*weightmap1f.ptr<float>(currv_nei)[w_idx]*cpp1d[w_idx].ptr<double>(tempFS+m)[tempFS+n];
									}
								}
							}
						}
					}else if(tempMinNumber==9){
						// update cross correlation
						for(int m=-tempFS;m<=tempFS;m++){
							for(int n=-tempFS;n<=tempFS;n++){
								if(i+m>=0&&i+m<m_Height&&j+n>=0&&j+n<m_Width){
									for(int w_idx=0;w_idx<2;w_idx++){
										c_ep1d[w_idx].ptr<double>(i+m)[j+n]+=a0[tempMinNumber]*weightmap1f.ptr<float>(currv)[w_idx]*cpp1d[w_idx].ptr<double>(tempFS+m)[tempFS+n];
									}
								}
							}
						}
					}
					BenefitPixelNumber++;
				}
			}
		}
		//		cout	<<	BenefitPixelNumber	<<	endl;
		if(BenefitPixelNumber==0){
			break;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	/// Change absorb to grayscale
	dst1b.create(src1b.size(),CV_8UC1);
	for(int i=0;i<m_Height;i++){
		for(int j=0;j<m_Width;j++){
			dst1b.ptr<uchar>(i)[j]=(1.-dst1f.ptr<float>(i)[j])*255.;
		}
	}

	return true;
}
