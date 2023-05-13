//
// Created by cs18 on 5/10/23.
//

#pragma once

#include <opencv2/opencv.hpp>
#include <iostream>
#include <cmath>
using namespace cv;
class DFT{
public:
    void My_DFT(Mat input_image, Mat& output_image, Mat& transform_image){
        //1.getOptimalDFTSize & expand size
        int m = getOptimalDFTSize(input_image.rows);
        int n = getOptimalDFTSize(input_image.cols);
        copyMakeBorder(input_image,input_image,0,
                       m-input_image.rows,0,n-input_image.cols,
                       BORDER_CONSTANT, Scalar::all(0));

        //2.To restore complex
        Mat planes[] = {Mat_<float>(input_image),
                        Mat::zeros(input_image.size(),CV_32F)};

        //3.merge planes into image
        merge(planes,2,transform_image);

        //4.DFT
        dft(transform_image,transform_image);

        //5.Get Magnitude
        split(transform_image,planes);
        magnitude(planes[0],planes[1],output_image);

        //6.log & regularize
        output_image += Scalar(1);
        log(output_image,output_image);
        normalize(output_image,output_image,0,1,NORM_MINMAX);

        //7.cut - make sure even
        output_image = output_image(Rect(0,0,output_image.cols & -2,
                                         output_image.rows & -2));

        //arrange
        int cx = output_image.cols / 2;
        int cy = output_image.rows / 2;
        Mat q0(output_image, Rect(0,0,cx,cy));
        Mat q1(output_image,Rect(cx,0,cx,cy));
        Mat q2(output_image, Rect(0,cy,cx,cy));
        Mat q3(output_image,Rect(cx,cy,cx,cy));

        //centralize
        Mat tmp;
        q0.copyTo(tmp);q3.copyTo(q0);tmp.copyTo(q3);
        q1.copyTo(tmp);q2.copyTo(q1);tmp.copyTo(q2);
    }

    void LPF_Filter(std::vector<double> & data, int cols, int rows, double sigma){
        float sigma22 = 2 * sigma * sigma;
        float temp = 0.0;
        int halfRows = cvRound(rows / 2);
        int halfCols = cvRound(cols / 2);
        char debug[12] = {0};
        for(int j = 1; j <= rows; ++j){
            for(int i = 1; i <= cols; ++i){
                temp = pow(j - halfRows, 2) + pow(i - halfCols, 2);
                data.push_back(exp(-temp / sigma22));
            }
        }
    }

    void genFreImg(std::string s){
        Mat image, image_gray, image_output, image_transform;
        image = imread(s);
        if(image.empty()){
            std::cout<<"READ ERROR!"<<std::endl;
            exit(-1);
        }

        cvtColor(image, image_gray, COLOR_BGR2GRAY);
        My_DFT(image_gray,image_output,image_transform);

        image_output *= 255;
        imwrite("output.ppm", image_output);
    }

    void genFreImg(std::string s1, std::string s2){
        Mat image, image_gray, image_output, image_transform;
        Mat imageA = imread(s1);
        if(imageA.empty()){
            std::cout<<"READ ERROR!"<<std::endl;
            exit(-1);
        }
        Mat imageB = imread(s2);
        if(imageB.empty()){
            std::cout<<"READ ERROR!"<<std::endl;
            exit(-1);
        }
        image = imageA - imageB;

        cvtColor(image, image_gray, COLOR_BGR2GRAY);

        My_DFT(image_gray,image_output,image_transform);

        image_output *= 255;
        imwrite("output.ppm", image_output);
    }

    void gaussFilter(Mat spec, Mat& res){
        std::vector<double>data;
        double sigma = sqrt(2 / acos(-1));
        LPF_Filter(data, spec.cols, spec.rows, sigma);
        spec.copyTo(res);

        for(size_t j = 0; j < spec.rows; ++j)
            for(size_t i = 0; i < spec.cols; ++i){
//                if(data.at(i + j * spec.cols) > 0.1) printf("%lf ",data.at(i + j * spec.cols)); if(i + 1 == spec.cols) printf("\n");
                *res.ptr<float>(j, i) *= data.at(i + j * spec.cols);
            }

        //arrange
        int cx = res.cols / 2;
        int cy = res.rows / 2;
        Mat q0(res, Rect(0,0,cx,cy));
        Mat q1(res,Rect(cx,0,cx,cy));
        Mat q2(res, Rect(0,cy,cx,cy));
        Mat q3(res,Rect(cx,cy,cx,cy));

        //decentralize
        Mat tmp;
        q0.copyTo(tmp);q3.copyTo(q0);tmp.copyTo(q3);
        q1.copyTo(tmp);q2.copyTo(q1);tmp.copyTo(q2);

    }

    void myFilter(std::string s){
        Mat image, image_output;
        image = imread(s);
        if(image.empty()){
            std::cout<<"READ ERROR!"<<std::endl;
            exit(-1);
        }

        double sigma = sqrt(2.0 / acos(-1));
        GaussianBlur(image, image_output, Size(3, 3), sigma, sigma);

        imwrite("filteredImg.ppm", image_output);
    }

    void tilewiseSpec(std::string s){
        Mat image, image_gray, image_output, image_transform;
        image = imread(s);
        if(image.empty()){
            std::cout<<"READ ERROR!"<<std::endl;
            exit(-1);
        }
        Mat GT = imread("sppp4D.ppm");
        image = image - GT;

        int tile = 32;
        copyMakeBorder(image,image,0,
                       tile-(image.rows%tile),0,tile-(image.cols%tile),
                       BORDER_CONSTANT, Scalar::all(0));

        cvtColor(image, image_gray, COLOR_BGR2GRAY);
        image_output = Mat::zeros(image_gray.size(), CV_32FC1);
        for(int i = 0; i < image.cols / tile; ++i)
            for(int j = 0; j < image.rows / tile; ++j){
                Mat qs(image_gray, Rect(i*tile, j*tile, tile, tile));
                Mat qd(image_output, Rect(i*tile, j*tile, tile, tile));
                My_DFT(qs,qd,image_transform);
            }
//        My_DFT(image_gray,image_output,image_transform);

        image_output *= 255;
        imwrite("output.ppm", image_output);
    }
};
