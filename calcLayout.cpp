#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <sys/time.h>
#include <Accelerate/Accelerate.h>

/*
extern "C" int ssyev_(char *jobz, char *uplo,
    int *n, float *a, int *lda, float *w,
    float *work, int *lwork, int *info);
*/
/*
extern "C" int dsyev_(char *jobz, char
        *uplo, int *n, double *a, int *lda, double *w,
        double *work, int *lwork, int *info);
*/

#define dimension 10000

typedef boost::numeric::ublas::vector<float> fvector;

using namespace std;

float * solver3D(float *, float *, float, float *);
float * solver2D(float *, float *, float, float *);

//Params about Graph
string graphName = "";
int N;
int M;
int ** D;
vector<int> * neighbor;
vector< pair<int, int> > edges;
vector<int> * edgelist;
vector<string> labels;
float * nodevalues; float * edgevalues;
float nodevalue_max, nodevalue_min;
float edgevalue_min, edgevalue_max;

//Graph Layouts
int dim;
float * lambdas;
float * P;
float * P_norms;
float * Layout3D;
float * E_3D;
float * E_3D_init;
float * Layout2D;
float * E_2D;
float * E_2D_init;

static float init3D[15];
static float init2D[6];

float scale = 1.0f;
float delta = 0.5;

static int pn_3d= 1; static int pn_2d= 1;
static int diff_3d= 1; static int diff_2d= 1;

static float memory_status = false;

double getETime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int strToInt(string &str){
    int t;
    stringstream ss;
    ss << str;
    ss >> t;
    return t;
}

string IntToString(int num){
    stringstream ss;
    ss << num;
    return ss.str();
}

void resetLayout3D(){
    scale = 1.0f;
    delta = 0.5;

    for(int i = 0; i < dim; i++){
        for(int j = 0; j < 3; j++){
            E_3D[i+j*dim] = E_3D_init[i+j*dim];
        }
    }

    delete[] Layout3D;
    Layout3D = new float[N*3];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 3, dim, scale, P, N, E_3D, dim, 0.0, Layout3D, N);

    {
        delete nodevalues;
        nodevalues = new float[N];
        for(int i = 0; i < N; i++){
            nodevalues[i] = 1.0f;
        }
        nodevalue_max = 5; nodevalue_min = 0;

        delete edgevalues;
        edgevalues = new float[M];
        for(int i = 0; i < M; i++){
            edgevalues[i] = 1.0f;
        }
        edgevalue_max = 1; edgevalue_min = 0;
    }
}

void resetLayout2D(){
    scale = 1.0f;
    delta = 0.5;

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < dim; j++){
            E_2D[j+i*dim] = E_2D_init[j+i*dim];
        }
    }

    delete[] Layout2D;
    Layout2D = new float[N*2];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 2, dim, scale, P, N, E_2D, dim, 0.0, Layout2D, N);

    {
        delete nodevalues;
        nodevalues = new float[N];
        for(int i = 0; i < N; i++){
            nodevalues[i] = 1.0f;
        }
        nodevalue_max = 5; nodevalue_min = 0;

        delete edgevalues;
        edgevalues = new float[M];
        for(int i = 0; i < M; i++){
            edgevalues[i] = 1.0f;
        }
        edgevalue_max = 1; edgevalue_min = 0;
    }
}

void calcmdsLayout(){
    delete[] lambdas;
    delete[] P_norms;
    delete[] P;

    float * D2 = new float[N*N];
    for(int i = 0; i < N; i++){
        D2[i+i*N] = 0.0f;
        for(int j = i+1; j < N; j++){
            float d2 = (float)(D[i][j]*D[i][j]);
            D2[i+j*N] = d2; D2[j+i*N] = d2;
        }
    }

    float * H = new float[N*N];
    for(int i = 0; i < N; i++){
        H[i+N*i] = 1.0f - (float)1.0/N;
        for(int j = i+1; j < N; j++){
            H[i+N*j] = -(float)1.0/N;
            H[j+N*i] = -(float)1.0/N;
        }
    }

    float * b = new float[N*N]; float * B = new float[N*N];
    float alpha = 1.0, beta = 0.0f;
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, N, N, alpha, H, N, D2, N, beta, b, N);
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, N, N, alpha, b, N, H, N, beta, B, N);

    double st2 = getETime();
    //Float
    int lwork, info;
    float *A = new float[N*N];
    float *w = new float[N];

    for(int i = 0; i < N; i++){
        A[i+i*N] = -0.5f*B[i+i*N];
        for(int j = i+1; j < N; j++){
            A[i+j*N] = -0.5f*B[i+j*N];
            A[j+i*N] = -0.5f*B[j+i*N];
        }
    }

    lwork = -1;
    float *work = new float[1];
    char jobz = 'V', uplo = 'U';
    ssyev_(&jobz, &uplo, &N, A, &N, w, work, &lwork, &info);
    lwork = (int)work[0];
    delete[]work;
    work = new float[max((int) 1, lwork)];
    ssyev_(&jobz, &uplo, &N, A, &N, w, work, &lwork, &info);

    //double
    /*
    int lwork, info;
    double *A = new double[N*N];
    double *w = new double[N];
    for(int i = 0; i < N; i++){
        A[i+i*N] = -0.5f*B[i+i*N];
        for(int j = i+1; j < N; j++){
            A[i+j*N] = -0.5f*B[i+j*N];
            A[j+i*N] = -0.5f*B[j+i*N];
        }
    }
    lwork = -1;
    double *work = new double[1];
    char jobz = 'V', uplo = 'U';
    dsyev_(&jobz, &uplo, &N, A, &N, w, work, &lwork, &info);
    lwork = (int)work[0];
    delete[]work;
    work = new double[max((int) 1, lwork)];
    dsyev_(&jobz, &uplo, &N, A, &N, w, work, &lwork, &info);
    */

    double en2 = getETime();
    cout << "eigensolver: " << en2 - st2 << endl;

    /*
    for(int i = 0; i < N; i++){
        cout << w[i] << endl;
    }*/

    dim = 0;
    float eps = 0.01;
    for(int i = N-1; i >= 0; i--){
        if(w[i] > eps) dim++;
        else break;
    }

    if(dim > dimension) dim = dimension;

    lambdas = new float[dim];
    for(int i = 0; i < dim; i++){
        lambdas[i] = w[N-1-i];
    }

    /*string filename = "lambdas.txt";
    ofstream ofs(filename.c_str());
    for(int i = 0; i < dim; i++){
        ofs << lambdas[i] << endl;
    }
    ofs.close();*/

    float * L = new float[N*dim];

    P = new float[N*dim];
    P_norms = new float[N];

    for(int i = 0; i < N; i++){
        for(int j = 0; j < dim; j++){
            if( i==j ) L[i+j*N] = (float)sqrt(w[N-1-j]);
            else L[i+j*N] = 0.0f;
            P[i+j*N] = 0.0f;
        }
    }

    float * U = new float[N*N];
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            U[i+j*N] = A[j+(N-1-i)*N];
        }
    }

    //High Dimensional Layout : P
    double st3 = getETime();
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
        N, dim, N, alpha, U, N, L, N, beta, P, N);
    double en3 = getETime();
    printf("projetion: %.6f\n", en3-st3);

    float pij = 0;
    for(int i = 0; i < N; i++){
        P_norms[i] = 0.0f;
        for(int j = 0; j < dim; j++){
            pij = P[i+j*N];
            P_norms[i] += pij*pij;
        }
    }

    delete work; delete w; delete A;
    delete D2; delete H; delete B; delete b;
    delete L; delete U;
}

void UpdateScale3D(float r){
    scale = r;
    delete[] Layout3D;
    Layout3D = new float[N*3];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 3, dim, scale, P, N, E_3D, dim, 0.0, Layout3D, N);
}

void UpdateScale2D(float r){
    scale = r;
    delete[] Layout2D;
    Layout2D = new float[N*2];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 2, dim,scale, P, N, E_2D, dim, 0.0, Layout2D, N);
}

void UpdateProjection3D(float r){
    delta = r;
    fvector e_3d[3];
    for(int i = 0; i < 3; i++){
        e_3d[i] = fvector(dim);
        for(int j = 0; j < dim; j++){
            if(j % 3 == i){
                e_3d[i][j] = pow((double)lambdas[j], (double)delta);
            }
            else e_3d[i][j] = 0.0f;
        }
        e_3d[i] = e_3d[i] / norm_2(e_3d[i]);
    }

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < dim; j++){
            E_3D[j+i*dim] = e_3d[i][j];
        }
    }

    delete[] Layout3D;
    Layout3D = new float[N*3];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 3, dim, scale, P, N, E_3D, dim, 0.0, Layout3D, N);
}

void UpdateProjection2D(float r){
    delta = r;
    fvector e_2d[2];
    for(int i = 0; i < 2; i++){
        e_2d[i] = fvector(dim);
        for(int j = 0; j < dim; j++){
            if(j % 2 == i) e_2d[i][j] = pow((double)lambdas[j], (double)delta);
            else e_2d[i][j] = 0.0f;
        }
        e_2d[i] = e_2d[i] / norm_2(e_2d[i]);
    }

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < dim; j++){
            E_2D[j+i*dim] = e_2d[i][j];
        }
    }

    delete[] Layout2D;
    Layout2D = new float[N*2];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 2, dim, scale, P, N, E_2D, dim, 0.0, Layout2D, N);
}

void UpdateDimension3D(float r){
    int _dim = (int)(r*dim);
    if(_dim < 3) _dim = 3;

    fvector e_3d[3];
    for(int i = 0; i < 3; i++){
        e_3d[i] = fvector(dim);
        for(int j = 0; j < dim; j++){
            if(j % 3 == i && j < _dim){
                e_3d[i][j] = pow((double)lambdas[j], (double)delta);
            }
            else e_3d[i][j] = 0.0f;
        }
        e_3d[i] = e_3d[i] / norm_2(e_3d[i]);
    }

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < dim; j++){
            E_3D[j+i*dim] = e_3d[i][j];
            E_3D_init[j+i*dim] = e_3d[i][j];
        }
    }

    delete[] Layout3D;
    Layout3D = new float[N*3];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 3, dim, scale, P, N, E_3D, dim, 0.0, Layout3D, N);
}

void UpdateDimension2D(float r){
    int _dim = (int)(r*dim);
    if(_dim < 2) _dim = 2;

    fvector e_2d[2];
    for(int i = 0; i < 2; i++){
        e_2d[i] = fvector(dim);
        for(int j = 0; j < dim; j++){
            if(j % 2 == i && j < _dim) e_2d[i][j] = pow((double)lambdas[j], (double)delta);
            else e_2d[i][j] = 0.0f;
        }
        e_2d[i] = e_2d[i] / norm_2(e_2d[i]);
    }

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < dim; j++){
            E_2D[j+i*dim] = e_2d[i][j];
        }
    }

    delete[] Layout2D;
    Layout2D = new float[N*2];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 2, dim, scale, P, N, E_2D, dim, 0.0, Layout2D, N);
}

void calc3DLayout(){
    fvector e[3];
    for(int i = 0; i < 3; i++){
        e[i] = fvector(dim);
        for(int j = 0; j < dim; j++){
            if(j % 3 == i){
                e[i][j] = sqrt(lambdas[j]);
            }
            else e[i][j] = 0.0f;
        }
        e[i] = e[i] / norm_2(e[i]);
    }

    E_3D = new float[dim*3];
    E_3D_init = new float[dim*3];

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < dim; j++){
            E_3D[j+i*dim] = e[i][j];
            E_3D_init[j+i*dim] = e[i][j];
        }
    }

    delete[] Layout3D;
    Layout3D = new float[N*3];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 3, dim, 1.0, P, N, E_3D, dim, 0.0, Layout3D, N);

    //Set Initial Value of resctriction
    init3D[0] = 0.9; init3D[1] = -0.1; init3D[2] = 0.1;
    init3D[3] = -0.1; init3D[4] = 0.9; init3D[5] = 0.1;
    init3D[6] = 0.1; init3D[7] = 0.1; init3D[8] = 0.9;
    init3D[9] = 0.1; init3D[10] = -0.1; init3D[11] = 0.1;
    init3D[12] = -0.1; init3D[13] = 0.1; init3D[14] = 0.1;
}

void calc2DLayout(){
    fvector e[2];
    for(int i = 0; i < 2; i++){
        e[i] = fvector(dim);
        for(int j = 0; j < dim; j++){
            if(j % 2 == i) e[i][j] = sqrt(lambdas[j]);
            else e[i][j] = 0.0f;
        }
        e[i] = e[i] / norm_2(e[i]);
    }

    E_2D = new float[dim*2];
    E_2D_init = new float[dim*2];

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < dim; j++){
            E_2D[j+i*dim] = e[i][j];
            E_2D_init[j+i*dim] = e[i][j];
        }
    }

    delete[] Layout2D;
    Layout2D = new float[N*2];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, 2, dim, scale, P, N, E_2D, dim, 0.0, Layout2D, N);

    //Set Initial Value of resctriction
    init2D[0] = 0.9; init2D[1] = -0.1;
    init2D[2] = 0.1; init2D[3] = 0.9;
    init2D[4] = 0.4; init2D[5] = 0.4;
}

void loadNodeAttrData(int n){
    string datadir = "../data/" + graphName + "/";
    bool isNormal = false;
    switch(n){
        case 1: datadir += graphName + "degree.txt"; break;
        case 2: datadir += graphName + "clcent.txt";  break;
        case 3: datadir += graphName + "bwcent.txt"; break;
        case 4: datadir += graphName + "clcoeff.txt"; break;
        case 5: datadir += graphName + "evcent.txt"; break;
        case 6: datadir += graphName + "pr.txt"; break;
        default: isNormal = true; break;
    }

    if(isNormal){
        for(int i = 0; i < N; i++){
            nodevalues[i] = 1.0f;
        }
        nodevalue_max = 5;
        nodevalue_min = 0;
    }
    else{
        ifstream ifs(datadir.c_str());
        if(!ifs.fail()){
            nodevalue_max = 0;
            nodevalue_min = 1000000000000;
            for(int i = 0; i < N; i++){
                ifs >> nodevalues[i];
                nodevalue_min = min(nodevalue_min, nodevalues[i]);
                nodevalue_max = max(nodevalue_max, nodevalues[i]);
            }
         }
        else{
            if(n == 1){
                nodevalue_max = 0;
                nodevalue_min = 1000000000000;
                for(int i = 0; i < N; i++){
                    nodevalues[i] = neighbor[i].size();
                    nodevalue_min = min(nodevalue_min, nodevalues[i]);
                    nodevalue_max = max(nodevalue_max, nodevalues[i]);
                }
            }
            else{
                cerr << "File not found\n";
            }
        }
        ifs.close();
    }
}

void loadEdgeAttrData(int n){
    string datadir = "../data/" + graphName + "/";
    bool isDefault = false;
    switch(n){
        case 1: datadir += graphName + "Simpson.txt"; break;
        case 2: datadir += graphName + "SimpsonEx.txt"; break;
        case 3: datadir += graphName + "ebw.txt"; break;
        case 4: datadir += graphName + "Weight.txt"; break;
        default: isDefault = true; break;
    }

    //cout << datadir << endl;

    if(isDefault){
        for(int i = 0; i < M; i++){
            edgevalues[i] = 1.0;
        }
        edgevalue_min = 0;
        edgevalue_max = 1;
    }
    else{
        ifstream ifs(datadir.c_str());
        if(!ifs.fail()){
            edgevalue_max= 0;
            edgevalue_min = 100000000;
            int from, to;
            for(int i = 0; i < M; i++){
                ifs >> from >> to >> edgevalues[i];
                edgevalue_min = min(edgevalue_min, edgevalues[i]);
                edgevalue_max = max(edgevalue_max, edgevalues[i]);
            }
            ifs.close();
        }
        else{
            cout << "File not found" << endl;
            //Calculate Attributes
            bool c = false;
            if(c){
                if(n == 1){
                    string filename = graphName + "Simpson.txt";
                    ofstream ofs(filename.c_str());
                    edgevalue_max= 0;
                    edgevalue_min = 100000000;
                    for(int i = 0; i < M; i++){
                        int from = edges[i].first, to = edges[i].second;
                        float inter = 1.0f;
                        for(int j = 0; j < N; j++){
                            if(D[from][j] == 1 && D[to][j] == 1){
                                inter += 1.0f;
                            }
                        }
                        float f = (float)min(neighbor[from].size(), neighbor[to].size());
                        edgevalues[i] = (float)inter / (float)f;
                        edgevalue_min = min(edgevalue_min, edgevalues[i]);
                        edgevalue_max = max(edgevalue_max, edgevalues[i]);
                        ofs << from << " " << to << " " << edgevalues[i] << endl;
                    }
                    ofs.close();
                }
                else if(n == 2){
                    string filename = graphName + "SimpsonEx.txt";
                    ofstream ofs(filename.c_str());
                    edgevalue_max= 0;
                    edgevalue_min = 100000000;
                    for(int i = 0; i < M; i++){
                        int from = edges[i].first, to = edges[i].second;
                        float count = 0.0f;
                        for(int j = 0; j < N; j++){
                            if(D[from][j] == 1 && D[to][j] == 1){
                                count += 1.0f;
                            }
                        }
                        int a = min(neighbor[from].size(), neighbor[to].size());
                        if(a >= 2){
                            edgevalues[i] = (float)count / (float)a;
                        }
                        else{
                            edgevalues[i] = 0;
                        }
                        edgevalue_min = min(edgevalue_min, edgevalues[i]);
                        edgevalue_max = max(edgevalue_max, edgevalues[i]);
                        ofs << from << " " << to << " " << edgevalues[i] << endl;
                    }
                    ofs.close();
                }
            }
            else{
                for(int i = 0; i < M; i++){
                    edgevalues[i] = 1.0;
                }
                edgevalue_min = 0;
                edgevalue_max = 1;
            }
        }
    }
}

void loadLabelData(){
    string labeldata = "../data/" + graphName + "/" + graphName + "labels.txt";
    ifstream lifs(labeldata.c_str());
    if(lifs.fail()){
        cout << "No Label Data" << endl;
        for(int i = 0; i < N; i++){
            labels.push_back(IntToString(i+1));
        }
    }
    else{
        string _str = "";

        for(int i = 0; i < N; i++){
            getline(lifs, _str);
            labels.push_back(_str);
        }
    }
    lifs.close();
}

//Set Default Node & Edge value
void setNodeEdgeValue(){
    delete[] nodevalues;
    nodevalues = new float[N];
    for(int i = 0; i < N; i++){
        nodevalues[i] = 1.0f;
    }
    nodevalue_max = 5; nodevalue_min = 0;

    delete[] edgevalues;
    edgevalues = new float[M];
    for(int i = 0; i < M; i++){
        edgevalues[i] = 1.0f;
    }
    edgevalue_max = 1; edgevalue_min = 0;
}

void loadMatrixData_t(const char * data){
    //Free Memory
    {
        if(memory_status){
            for(int i = 0; i < N; i++){
                delete[] D[i];
            }
            delete[] D;
            memory_status = false;
        }
        for(int i = 0; i < N; i++){
            vector<int>().swap(neighbor[i]);
            vector<int>().swap(edgelist[i]);
        }
        vector< pair<int, int> >().swap(edges);
        vector<string>().swap(labels);
    }

    {
        double load_s = getETime();
        ifstream ifs(data);
        if(ifs.fail()){
            cerr << "File not found\n";
            exit(0);
        }

        ifs >> N;
        D = new int*[N];
        memory_status = true;
        neighbor = new vector<int>[N];
        edgelist = new vector<int>[N];

        for(int i = 0; i < N; i++){
            D[i] = new int[N];
            for(int j = 0; j < N; j++){
                ifs >> D[i][j];
            }
        }

        int diameter = 0;

        for(int i = 0; i < N; i++){
            for(int j = i+1; j < N; j++){
                diameter = max(diameter, D[i][j]);
                if(D[i][j] == 1){
                    neighbor[i].push_back(j);
                    neighbor[j].push_back(i);
                    edges.push_back(make_pair(i,j));
                }
            }
        }

        M = edges.size();
        for(int i = 0; i < M; i++){
            int from = edges[i].first, to = edges[i].second;
            //cout << from << " " << to << endl;
            edgelist[from].push_back(i);
            edgelist[to].push_back(i);
        }

        double load_e = getETime();
        cout << "Load Time: " << load_e - load_s << endl;
        double mds_s = getETime();

        calcmdsLayout();
        //double mds_e = getETime();
        //cout << "MDS Time: " << mds_e - mds_s << endl;
        //double pro_s = getETime();
        calc3DLayout();
        double pro_e = getETime();
        //cout << "Projection Time: " << pro_e - pro_s << endl;
        cout << "Layout Calc Time: " << pro_e - mds_s << endl;
        calc2DLayout();

        cout << N << " " << M << " " << diameter << " " << dim << endl;
        ifs.close();
    }

    //Load Label Data
    loadLabelData();
    //Set Default Node & Edge value
    setNodeEdgeValue();

    pn_3d= 1; pn_2d= 1; diff_3d= 1; diff_2d= 1;
    scale = 1.0f; delta = 0.5;
}

void loadMatrixData_b(const char * data){
    //Free Memory
    {
        if(memory_status){
            for(int i = 0; i < N; i++){
                delete[] D[i];
            }
            delete[] D;
            memory_status = false;
        }
        for(int i = 0; i < N; i++){
            vector<int>().swap(neighbor[i]);
            vector<int>().swap(edgelist[i]);
        }
        vector< pair<int, int> >().swap(edges);
        vector<string>().swap(labels);
    }

    {
        //double load_s = getETime();

        ifstream ifs(data, ios::out|ios::binary);
        if(ifs.fail()){
            cerr << "File not found\n";
            exit(0);
        }

        int _N;
        ifs.read((char *)&_N, sizeof(int));
        N = _N;
        int * buffer = (int *)malloc((long)sizeof(int)*(N*N));
        ifs.read((char *)buffer, (long)sizeof(int)*N*N);
        ifs.close();

        D = new int*[N];
        memory_status = true;
        neighbor = new vector<int>[N];
        edgelist = new vector<int>[N];

        for(int i = 0; i < N; i++){
            D[i] = new int[N];
            for(int j = 0; j < N; j++){
                D[i][j] = buffer[i+j*N];
            }
        }
        free(buffer);

        int diameter = 0;
        for(int i = 0; i < N; i++){
            for(int j = i+1; j < N; j++){
                diameter = max(diameter, D[i][j]);
                if(D[i][j] == 1){
                    neighbor[i].push_back(j);
                    neighbor[j].push_back(i);
                    edges.push_back(make_pair(i,j));
                }
            }
        }
        M = edges.size();

        for(int i = 0; i < M; i++){
            int from = edges[i].first, to = edges[i].second;
            edgelist[from].push_back(i);
            edgelist[to].push_back(i);
        }

        //double load_e = getETime();
        //cout << "Load Time: " << load_e - load_s << endl;

        double mds_s = getETime();
        calcmdsLayout();
        //double mds_e = getETime();
        //cout << "MDS Time: " << mds_e - mds_s << endl;
        //double pro_s = getETime();
        calc3DLayout();
        double pro_e = getETime();
        //cout << "Projection Time: " << pro_e - pro_s << endl;
        cout << "Layout Calc Time: " << pro_e - mds_s << endl;
        calc2DLayout();

        cout << N << " " << M << " " << diameter << " " << dim << endl;
        ifs.close();
    }

    //Load Label Data
    loadLabelData();
    //Set Default Node & Edge value
    setNodeEdgeValue();
    pn_3d= 1; pn_2d= 1; diff_3d= 1; diff_2d= 1;
    scale = 1.0f; delta = 0.5;
}

void loadLayoutData_t(const char * data){
    //Free Memory
    {
        if(memory_status){
            for(int i = 0; i < N; i++){
                delete[] D[i];
            }
            delete[] D;
            memory_status = false;
        }
        for(int i = 0; i < N; i++){
            vector<int>().swap(neighbor[i]);
            vector<int>().swap(edgelist[i]);
        }
        vector< pair<int, int> >().swap(edges);
        vector<string>().swap(labels);
        delete[] lambdas;
        delete[] P_norms;
        delete[] P;
    }

    {
        double load_s = getETime();

        ifstream ifs(data);
        string str;

        if(ifs.fail()){
            cerr << "File not found\n";
            exit(0);
        }

        ifs >> str >> N >> M;

        neighbor = new vector<int>[N];
        edgelist = new vector<int>[N];

        //load Edges??? and binary
        for(int i = 0; i < M; i++){
            int _from, _to;
            ifs >> _from >> _to;
            neighbor[_from].push_back(_to); neighbor[_to].push_back(_from);
            edges.push_back(make_pair(_from,_to));
            edgelist[_from].push_back(i); edgelist[_to].push_back(i);
        }

        ifs >> dim;

        lambdas = new float[dim];
        for(int i = 0; i < dim; i++){
            ifs >> lambdas[i];
        }

        P = new float[N*dim];
        P_norms = new float[N];

        float pij = 0;
        for(int i = 0; i < N; i++){
            P_norms[i] = 0;
            for(int j = 0; j < dim; j++){
                ifs >> pij;
                P[i+j*N] = pij;
                P_norms[i] += pij*pij;
            }
        }

        double load_e = getETime();
        cout << "Load Time: " << load_e - load_s << endl;

        double pro_s = getETime();
        calc3DLayout();
        calc2DLayout();
        double pro_e = getETime();
        cout << "Projection Time: " << pro_e - pro_s << endl;

        cout << N << " " << M << " " << dim << endl;
        ifs.close();
    }

    //Load Label Data
    loadLabelData();
    //Set Default Node & Edge value
    setNodeEdgeValue();

    pn_3d= 1; pn_2d= 1; diff_3d= 1; diff_2d= 1;
    scale = 1.0f; delta = 0.5;
}

void loadLayoutData_b(const char * data){
    //Free Memory
    {
        if(memory_status){
            for(int i = 0; i < N; i++){
                delete[] D[i];
            }
            delete[] D;
            memory_status = false;
        }
        for(int i = 0; i < N; i++){
            vector<int>().swap(neighbor[i]);
            vector<int>().swap(edgelist[i]);
        }
        vector< pair<int, int> >().swap(edges);
        vector<string>().swap(labels);
        delete[] lambdas;
        delete[] P_norms;
        delete[] P;
    }

    {
        double load_s = getETime();
        ifstream layoutData(data, ios::out|ios::binary);

        if(layoutData.fail()){
            cerr << "File not found\n";
            exit(0);
        }

        int _N = 0, _dim = 0;
        layoutData.read((char *)&_N, sizeof(int));
        layoutData.read((char *)&_dim, sizeof(int));
        N = _N; dim = _dim;

        float * _lambdas = new float[dim];
        layoutData.read((char *)_lambdas, (long)sizeof(float)*(dim));
        lambdas = new float[dim];
        for(int i = 0; i < dim; i++){
            lambdas[i] = _lambdas[i];
        }
        free(_lambdas);

        P = new float[(long)(N)*(long)(dim)];
        layoutData.read((char *)P, (long)sizeof(float)*((long)(N)*(long)(dim)));

        P_norms = new float[N];
        float pij = 0;
        for(int i = 0; i < N; i++){
            P_norms[i] = 0;
            for(int j = 0; j < dim; j++){
                pij = P[i+j*N];
                P_norms[i] += pij*pij;
            }
        }
        double load_e = getETime();
        cout << "Load Time: " << load_e - load_s << endl;

        double pro_s = getETime();
        calc3DLayout();
        calc2DLayout();
        double pro_e = getETime();
        cout << "Projection Time: " << pro_e - pro_s << endl;

        layoutData.close();
    }

    //Load Edges
    {
        string edgefile = "../data/" + graphName + "/" + graphName + "edges.bin";
        ifstream edgeData(edgefile.c_str(), ios::out|ios::binary);
        if(edgeData.fail()){
            cerr << "File not found\n";
            exit(0);
        }
        int _M = 0;
        edgeData.read((char *)&_M, sizeof(int));
        cout << _M << endl;
        M = _M;
        //D = new int*[N];
        //for(int i = 0; i < N; i++){
        //    D[i] = new int[N];
        //}
        int * _edges = new int[2*M];
        edgeData.read((char *)_edges, (long)sizeof(int)*(2*M));
        neighbor = new vector<int>[N];
        edgelist = new vector<int>[N];
        for(int i = 0; i < M; i++){
            int _from = _edges[2*i]; int _to = _edges[2*i+1];
            //D[_from][_to] = 1; D[_to][_from] = 1;
            neighbor[_from].push_back(_to);
            neighbor[_to].push_back(_from);
            edges.push_back(make_pair(_from,_to));
            edgelist[_from].push_back(i);
            edgelist[_to].push_back(i);
        }
        free(_edges);
        edgeData.close();
    }
    //Load Label Data
    loadLabelData();
    //Set Default Node & Edge value
    setNodeEdgeValue();
    cout << N << " " << M << " " << dim << endl;

    pn_3d= 1; pn_2d= 1; diff_3d= 1; diff_2d= 1;
    scale = 1.0f; delta = 0.5;
}

int getNew3DLayout(int id,
    float pre_x, float pre_y, float pre_z,
    float new_x, float new_y, float new_z){

    float _pre[3]; float _new[3];
    _pre[0] = pre_x; _pre[1] = pre_y; _pre[2] = pre_z;
    _new[0] = new_x; _new[1] = new_y; _new[2] = new_z;

    fvector p(dim);
    float p_norm = P_norms[id]*scale, new_norm = 0, pre_norm = 0;

    for(int i = 0; i < dim; i++){
        p(i) = scale*P[id+i*N];
    }

    for(int i = 0; i < 3; i++){
        new_norm += _new[i]*_new[i];
        pre_norm += _pre[i]*_pre[i];
    }

    if(new_norm < p_norm*0.95f && pre_norm < p_norm*0.95f){
        /*
           float t = 0.001;
           if(new_norm > (1-t)*p_norm){
                for(int i = 0; i < 3; i++){
                _new[i] *= (1-t)*p_norm/new_norm;
                }
           }
        */
           fvector e0(dim), f0(dim), e1(dim), e2(dim), e3(dim);

           for(int i = 0; i < dim; i++){
            e1(i) = E_3D[i+0*dim];
            e2(i) = E_3D[i+1*dim];
            e3(i) = E_3D[i+2*dim];
        }

        f0 = p - _pre[0]*e1 - _pre[1]*e2 - _pre[2]*e3;
        float norm_f0 = norm_2(f0);
        e0 = f0 / norm_f0;

        //cout << inner_prod(e0,e1) << " " << inner_prod(e0,e2) << " " << inner_prod(e0,e3) <<  endl;
        //cout << norm_f0 << endl;
        //double _s = getETime();
        float * res = solver3D(_pre, _new, norm_f0, init3D);
        //double _e = getETime();
        //cout << "solver:" << _e - _s << endl;

        if(res == NULL){
            init3D[0] = 0.9; init3D[1] = -0.1; init3D[2] = 0.1;
            init3D[3] = -0.1; init3D[4] = 0.9; init3D[5] = 0.1;
            init3D[6] = 0.1; init3D[7] = 0.1; init3D[8] = 0.9;
            init3D[9] = 0.4; init3D[10] = -0.4; init3D[11] = 0.4;
            init3D[12] = -0.4; init3D[13] = 0.4; init3D[14] = 0.4;
            cout << "NULL" << endl;
            return 0;
        }

        fvector e1_new(dim), e2_new(dim), e3_new(dim), r1(dim), r2(dim);

        float a10 = (new_x - _pre[0]*res[0] - _pre[1]*res[1] - _pre[2]*res[2])/norm_f0;
        float a20 = (new_y - _pre[0]*res[3] - _pre[1]*res[4] - _pre[2]*res[5])/norm_f0;
        float a30 = (new_z - _pre[0]*res[6] - _pre[1]*res[7] - _pre[2]*res[8])/norm_f0;

        e1_new = a10*e0 + res[0]*e1 + res[1]*e2 + res[2]*e3;
        e2_new = a20*e0 + res[3]*e1 + res[4]*e2 + res[5]*e3;
        e3_new = a30*e0 + res[6]*e1 + res[7]*e2 + res[8]*e3;

        r1 = res[9]*e1 + res[10]*e2 + res[11]*e3;
        r2 = res[12]*e1 + res[13]*e2 + res[14]*e3;

        /*
        float x = inner_prod(p,e1_new), y = inner_prod(p,e2_new), z = inner_prod(p,e3_new);
        float n_1 = norm_2(e1_new), n_2 = norm_2(e2_new), n_3 = norm_2(e3_new),
        n_r1 = norm_2(r1), n_r2 = norm_2(r2);
        float er1 = inner_prod(e1,r1) - inner_prod(e1_new,r1);
        float er2 = inner_prod(e1,r2) - inner_prod(e1_new,r2);
        float er3 = inner_prod(e2,r1) - inner_prod(e2_new,r1);
        float er4 = inner_prod(e2,r2) - inner_prod(e2_new,r2);
        float er5 = inner_prod(e3,r1) - inner_prod(e3_new,r1);
        float er6 = inner_prod(e3,r2) - inner_prod(e3_new,r2);
        float e1_e2 = inner_prod(e1_new,e2_new),
        e2_e3 = inner_prod(e2_new,e3_new),
        e1_e3 = inner_prod(e1_new,e3_new),
        r1_r2 = inner_prod(r1,r2);
        */
        /*
        float errors1 = 0.0f;
        errors1 += (x - _new[0])*(x - _new[0]);
        errors1 += (y - _new[1])*(y - _new[1]);
        errors1 += (z - _new[2])*(z - _new[2]);
        float errors2 = 0.0f;
        errors2 += (n_1 - 1)*(n_1 - 1);
        errors2 += (n_2 - 1)*(n_2 - 1);
        errors2 += (n_3 - 1)*(n_3 - 1);
        errors2 += (n_r2 - 1 )*(n_r2 - 1);
        float errors3 = 0.0f;
        errors3 += (e1_e2)*(e1_e2);
        errors3 += (e2_e3)*(e2_e3);
        errors3 += (e1_e3)*(e1_e3);
        errors3 += (r1_r2)*(r1_r2);
        float errors4 = 0.0f;
        errors4 += er1*er1;
        errors4 += er2*er2;
        errors4 += er3*er3;
        errors4 += er4*er4;
        errors4 += er5*er5;
        errors4 += er6*er6;
        cout << errors1 << endl;
        cout << errors2 << endl;
        cout << errors3 << endl;
        cout << errors4 << endl;
        cout << sqrt(errors1+errors2+errors3+errors4) << endl;
        cout << "== debug ==" << endl;
        cout << x - _new[0] << endl;
        cout << y - _new[1] << endl;
        cout << z - _new[2] << endl;
        cout << n_1 - 1 << endl;
        cout << n_2 - 1<< endl;
        cout << n_3 - 1<< endl;
        cout << n_r1 - 1 << endl;
        cout << n_r2 - 1 << endl;
        cout << e1_e2 << endl;
        cout << e2_e3 << endl;
        cout << e1_e3 << endl;
        cout << r1_r2 << endl;
        */

        float error = inner_prod(e1,e1_new) + inner_prod(e2,e2_new) + inner_prod(e3,e3_new);
        if(error < 2.8) return 0;

        /*Gram-Schmidt orthonormalization*/
        e1_new = e1_new / norm_2(e1_new);
        e2_new = e2_new - inner_prod(e1_new,e2_new)*e1_new;
        e2_new = e2_new / norm_2(e2_new);
        e3_new = e3_new - inner_prod(e1_new,e3_new)*e1_new - inner_prod(e2_new,e3_new)*e2_new;;
        e3_new = e3_new / norm_2(e3_new);

        for(int i = 0; i < dim; i++){
            E_3D[i+0*dim] = e1_new(i);
            E_3D[i+1*dim] = e2_new(i);
            E_3D[i+2*dim] = e3_new(i);
        }

        delete[] Layout3D;
        Layout3D = new float[N*3];
        double a = getETime();
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, 3, dim, scale, P, N, E_3D, dim, 0.0, Layout3D, N);
        double b = getETime();
        //cout << b - a << endl;
        return 1;
    }
    else{
        return 0;
    }
}

int getNew2DLayout(int id,
    float pre_x, float pre_y,
    float new_x, float new_y){

    float _pre[2]; float _new[2];
    _pre[0] = pre_x; _pre[1] = pre_y;
    _new[0] = new_x; _new[1] = new_y;

    fvector p(dim);
    float p_norm = P_norms[id]*scale, new_norm = 0, pre_norm = 0;

    for(int i = 0; i < dim; i++){
        p(i) = scale*P[id+i*N];
    }

    for(int i = 0; i < 2; i++){
        new_norm += _new[i]*_new[i];
        pre_norm += _pre[i]*_pre[i];
    }

    if(new_norm < p_norm*0.95f && pre_norm < p_norm*0.95f){
        /*
           float t = 0.001;

           if(new_norm > (1-t)*p_norm){
           for(int i = 0; i < 3; i++){
           _new[i] *= (1-t)*p_norm/new_norm;
           }
           }
           */
           fvector e0(dim), f0(dim), e1(dim), e2(dim);

           for(int i = 0; i < dim; i++){
            e1(i) = E_2D[i+0*dim];
            e2(i) = E_2D[i+1*dim];
        }

        f0 = p - _pre[0]*e1 - _pre[1]*e2;
        float norm_f0 = norm_2(f0);
        e0 = f0 / norm_f0;
        //cout << inner_prod(e0,e1) << " " << inner_prod(e0,e2) << " " << inner_prod(e0,e3) <<  endl;
        //cout << norm_f0 << endl;

        float * res = solver2D(_pre, _new, norm_f0, init2D);

        if(res == NULL){
            init2D[0] = 0.9; init2D[1] = -0.1;
            init2D[2] = 0.1; init2D[3] = 0.9;
            init2D[4] = 0.4; init2D[5] = 0.4;
            cout << "NULL" << endl;
            return 0;
        }

        fvector e1_new(dim), e2_new(dim), r(dim);

        float a10 = (new_x - _pre[0]*res[0] - _pre[1]*res[1])/norm_f0;
        float a20 = (new_y - _pre[0]*res[2] - _pre[1]*res[3])/norm_f0;

        e1_new = a10*e0 + res[0]*e1 + res[1]*e2;
        e2_new = a20*e0 + res[2]*e1 + res[3]*e2;

        //init
        /*
           for(int i = 0; i < 15; i++){
           init[i] = res[i];
           }
           */
           r = res[4]*e1 + res[5]*e2;

        /*
           float x = inner_prod(p,e1_new), y = inner_prod(p,e2_new);
           float n_1 = norm_2(e1_new), n_2 = norm_2(e2_new), n_r = norm_2(r);
           float er1 = inner_prod(e1,r) - inner_prod(e1_new,r);
           float er2 = inner_prod(e2,r) - inner_prod(e2_new,r);
           float e1_e2 = inner_prod(e1_new,e2_new);
         */
        /*
        cout << "== debug ==" << endl;
        cout << x - _new[0] << endl;
        cout << y - _new[1] << endl;
        cout << n_1 - 1 << endl;
        cout << n_2 - 1<< endl;
        cout << n_r - 1 << endl;
        */

        float error = inner_prod(e1,e1_new) + inner_prod(e2,e2_new);
        if(error < 1.6) return 0;

        //Gram-Schmidt orthonormalization
        e1_new = e1_new / norm_2(e1_new);
        e2_new = e2_new - inner_prod(e1_new,e2_new)*e1_new;
        e2_new = e2_new / norm_2(e2_new);

        for(int i = 0; i < dim; i++){
            E_2D[i+0*dim] = e1_new(i);
            E_2D[i+1*dim] = e2_new(i);
        }

        delete[] Layout2D;
        Layout2D = new float[N*2];
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, 2, dim,
            scale, P, N, E_2D, dim, 0.0, Layout2D, N);

        return 1;
    }
    else{
        return 0;
    }
}

void printLayout3D(){
    string filename = graphName + "Layout3D" + IntToString(pn_3d) + ".txt";
    ofstream ofs(filename.c_str());
    for(int i = 0; i < N; i++){
        ofs << Layout3D[i+0*N] << " " << Layout3D[i+1*N] << " " << Layout3D[i+2*N] << endl;
    }
    ofs.close();
    pn_3d++;
}

void printLayout2D(){
    string filename = graphName + "Layout2D" + IntToString(pn_2d) + ".txt";
    ofstream ofs(filename.c_str());
    for(int i = 0; i < N; i++){
        ofs << Layout2D[i+0*N] << " " << Layout2D[i+1*N] << endl;;
    }
    ofs.close();
    pn_2d++;
}

void calcStress3D(){
    string filename = graphName + "Diff3d" + IntToString(diff_3d) + ".txt";
    ofstream ofs(filename.c_str());
    float s = 0;
    float d = 0;
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            float x = (Layout3D[i+0*N] - Layout3D[j+0*N])/scale;
            float y = (Layout3D[i+1*N] - Layout3D[j+1*N])/scale;
            float z = (Layout3D[i+2*N] - Layout3D[j+2*N])/scale;
            float dij = sqrt( x*x + y*y + z*z );
            ofs << D[i][j] << " " << dij << endl;
            s += (dij - D[i][j])*(dij - D[i][j]);
            d += D[i][j]*D[i][j];
        }
    }
    ofs.close();
    diff_3d++;
    cout << "delta " << delta << endl;
    cout << "Stress Value " << s << endl;
    cout << "Normalized Stress Value " << sqrt(s / d) << endl;
}

void calcStress2D(){
    string filename = graphName + "Diff2d" + IntToString(diff_2d) + ".txt";
    ofstream ofs(filename.c_str());
    float s = 0;
    float d = 0;
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            float x = (Layout2D[i+0*N] - Layout2D[j+0*N])/scale;
            float y = (Layout2D[i+1*N] - Layout2D[j+1*N])/scale;
            float dij = sqrt( x*x + y*y );
            ofs << D[i][j] << " " << dij << endl;
            s += (dij - D[i][j])*(dij - D[i][j]);
            d += D[i][j]*D[i][j];
        }
    }
    ofs.close();
    diff_2d++;
    cout << delta << " " << sqrt(s / d) << endl;
}