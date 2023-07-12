import scipy.io as scio
import numpy as np
from sklearn.cross_decomposition import PLSRegression
import pycaltransfer.caltransfer as caltransfer
import dipals as ml
import functions as algo
import matplotlib.pyplot as plt

def predict_SST(Xt,F,a,BETA): 
    Xt_hat =   Xt.dot(F)+a
    col = np.ones(Xt.shape[0])
    Xt_hat = np.insert(Xt_hat,0,col,axis=1)
    Yhat = Xt_hat.dot(BETA) 
    return Yhat

def predict_DOP(Xs,Ys,Xt,E,LV): 
    Xs_e = Xs.dot(E)
    x_mean_e = np.mean(Xs_e, axis = 0)
    x_mean_e.shape = (1,Xs_e.shape[1])
    y_mean = np.mean(Ys, axis = 0)
    # PLSR after EPO
    pls2 = PLSRegression(n_components=LV,scale=False)
    pls2.fit(Xs_e, y)
    B_s = pls2.coef_
    beta_s = y_mean - (x_mean_e.dot(B_s)) # Model after TOP y = X_s B_s + beta_s
    Yhat = Xt.dot(B_s)+ beta_s
    return Yhat    
def predict_uDOP(Xs,Ys,Xt,E,LV): 
    Xs_e = Xs.dot(E) 
    x_mean_e = np.mean(Xs_e, axis = 0)
    x_mean_e.shape = (1,Xs_e.shape[1])
    y_mean = np.mean(Ys, axis = 0)
    pls2 = PLSRegression(n_components=LV,scale=False)
    pls2.fit(Xs_e, Ys)
    B = pls2.coef_
    beta = y_mean - (x_mean_e.dot(B)) # Model after TOP y = X_s B_s + beta_s
    Yhat = (Xt.dot(B) + beta) # predicted values in target domain
    return Yhat
#def evaluate_data():
#ReadData 
data = scio.loadmat('./dataset.mat')
TrainXs = data['TrainXs']
TrainYs = data['TrainYs']
TrainXt = data['TrainXt']
TrainYt = data['TrainYt']
TestXs = data['TestXs']
TestXt = data['TestXt']
TestYs = data['TestYs']
TestYt = data['TestYt']
indices = data['indices']
fold = data['fold'][0][0]
LV = data['LV'][0][0]
rhoArray = data['rhoArray'][0]
espArray = data['espArray'][0]
compArray = data['compArray'][0]
BETA = data['BETA']
tasknum = data['tasknum'][0]
task = ''
for i in range(0,tasknum.shape[0]):
    task += chr(tasknum[i])

parameter = np.array([rhoArray.size,espArray.size,compArray.size])
max_p = np.max(parameter)
RMSECV_DIPLS = np.zeros((rhoArray.size, fold), dtype=float)
RMSECV_SST = np.zeros((compArray.size, fold), dtype=float)
RMSECV_DOP = np.zeros((espArray.size, fold), dtype=float)

p = TrainYs.shape[1]
RMSEP_uDOP = np.zeros((compArray.size,p), dtype=float)
RMSEP_DOP  = np.zeros((1,p), dtype=float)
RMSEP_DIPLS  = np.zeros((1,p), dtype=float)
RMSEP_SST  = np.zeros((1,p), dtype=float)
RDP_uDOP = np.zeros((compArray.size,p), dtype=float)
RDP_DOP  = np.zeros((1,p), dtype=float)
RDP_DIPLS  = np.zeros((1,p), dtype=float)
RDP_SST  = np.zeros((1,p), dtype=float)


for i in range(0, max_p):
    for j in range(0, fold):
        CvIdx = np.where(indices != (j+1))
        CvIdx = CvIdx[0]
        CvTrainXs = TrainXs[CvIdx, :]
        CvTrainYs = TrainYs[CvIdx, :]
        CvTrainXt = TrainXt[CvIdx, :]
        CvTrainYt = TrainYt[CvIdx, :]
        CvIdx = np.where(indices == (j+1))
        CvIdx = CvIdx[0]
        CvTestXs = TrainXs[CvIdx, :]
        CvTestYs = TrainYs[CvIdx, :]
        CvTestXt = TrainXt[CvIdx, :]
        CvTestYt = TrainYt[CvIdx, :]
        #Evaluated by SST
        if i<compArray.size:
            comp = compArray[i]
            F, a = caltransfer.sst(CvTrainXs, CvTrainXt, ncomp = comp)
            Y_SST = predict_SST(CvTestXt,F,a,BETA)
            for k in range(0, p):
                RMSECV_SST[i,j] =  RMSECV_SST[i,j]+algo.rmse(Y_SST[:,k],CvTestYt[:,k])
             
        for k in range(0, p):
            y = CvTrainYs[:, k]
            y = y.reshape(y.shape[0], 1)
            yt = CvTrainYt[:, k]
            yt = yt.reshape(yt.shape[0], 1)
            ytest = CvTestYt[:, k]
            ytest = ytest.reshape(ytest.shape[0], 1)
             #Evaluated by di-pls  
            if i<rhoArray.size:
                rho = rhoArray[i]
                m = ml.model(CvTrainXs, y, CvTrainXs, CvTrainXt, LV)
                m.fit(int(rho))
                yhat_dipls, err = m.predict(CvTestXt, y_test=ytest)
                RMSECV_DIPLS[i, j] = RMSECV_DIPLS[i, j] + algo.rmse(yhat_dipls, ytest)  
            #Evaluated by DOP
            if i<espArray.size:
                E = caltransfer.dop(CvTrainXs, y, CvTrainXt, yt, rho = espArray[i], dop_ncomp = LV)
                Yhat = predict_DOP(CvTrainXs,y,CvTestXt,E,LV)
                RMSECV_DOP[i,j] =  RMSECV_DOP[i,j]+algo.rmse(Yhat,ytest)
    #Evaluated by uDOP
    if i<compArray.size:
        E = caltransfer.udop(TrainXs,  TrainXt, udop_ncomp = compArray[i])
        Yhat = predict_uDOP(TrainXs,TrainYs,TestXt,E,LV)
        for k in range(0,p):
            RMSEP_uDOP[i,k] = algo.rmse(Yhat[:,k],TestYt[:,k])
            RDP_uDOP[i,k] = np.std(TestYt[:,k])/ RMSEP_uDOP[i,k]
#calculate RMSEP  
#DIPLS optimal parameters
avg = np.average(RMSECV_DIPLS, axis=1)
i = np.argmin(avg)
rho1 =  rhoArray[i]
#DOP optimal parameters
avg = np.average(RMSECV_DOP, axis=1)
i = np.argmin(avg)
rho2 =  espArray[i]
#SST optimal parameters
avg = np.average(RMSECV_SST, axis=1)
i = np.argmin(avg)
comp =  compArray[i]
DIPLS_rho = rho1
DOP_rho = rho2
SST_comp = comp
for k in range(0, p):
    y = TrainYs[:, k]
    y = y.reshape(y.shape[0], 1)
    yt = TrainYt[:, k]
    yt = yt.reshape(yt.shape[0], 1)
    ytest = TestYt[:, k]
    ytest = ytest.reshape(ytest.shape[0], 1)
    #DIPLS
    #TrainXSel = 
    m = ml.model(TrainXs, y, TrainXs[1:int(TrainXs.shape[0]/2),:], TrainXt[int(TrainXs.shape[0]/2):TrainXs.shape[0],:], LV)
    m.fit(int(rho1))
    yhat_dipls, err = m.predict(TestXt, y_test=ytest)
    error = algo.rmse(yhat_dipls, ytest)
    RMSEP_DIPLS[0,k] = error
    RDP_DIPLS[0,k] = np.std(ytest)/error
    #DOP
    E = caltransfer.dop(TrainXs,y,TrainXt, yt, rho = rho2, dop_ncomp = LV)
    Yhat = predict_DOP(TrainXs,y,TestXt,E,LV)
    RMSEP_DOP[0,k] = algo.rmse(Yhat,ytest)
    RDP_DOP[0,k] =  np.std(ytest)/RMSEP_DOP[0,k]
#SST
F, a = caltransfer.sst(TrainXs, TrainXt, ncomp = comp)
Y_SST = predict_SST(TestXt,F,a,BETA)
for k in range(0,p):
    RMSEP_SST[0,k] =  algo.rmse(Y_SST[:,k],TestYt[:,k])
    RDP_SST[0,k] =  np.std(TestYt[:,k])/RMSEP_SST[0,k]
#RMSEP_SST /= p
#RMSEP_DOP /= p
#RMSEP_DIPLS /= p
#RMSEP_uDOP/= p

scio.savemat('./'+task+'_result.mat',{'RMSEP_SST':RMSEP_SST,\
    'RMSEP_DOP':RMSEP_DOP,'RMSEP_DIPLS':RMSEP_DIPLS,'RMSEP_uDOP':RMSEP_uDOP,\
    'RDP_SST':RDP_SST,'RDP_DOP':RDP_DOP,'RDP_DIPLS':RDP_DIPLS,\
        'RDP_uDOP':RDP_uDOP,'DIPLS_rho':DIPLS_rho,'DOP_rho':DOP_rho,'SST_comp':SST_comp})
 