import sys
import numpy as np
import os
import scipy.io as io
import torch
import torch.nn.functional as F
from time import time
from scipy.io import savemat
from PIL import Image


def tiff3Dread(path):
    img = Image.open(path)
    images = []
    for i in range(img.n_frames):
        img.seek(i)
        images.append(np.array(img))
    return np.array(images)


class parameter:
    def __init__(self):
        # Phantom
        self.PhantomSize = []


class Modle_A_G(torch.nn.Module):

    def __init__(self, param, N, M):
        super().__init__()
        Num_Ray = [N, M]
        # Initialize_Refraction
        ShiftMap = Initialize_ShiftMap_ept(Num_Ray)
        ShiftMap = ShiftMap[np.newaxis, :]
        ShiftMap = ShiftMap.astype('float32')
        ShiftMap = torch.tensor(ShiftMap)
        ShiftMap_gpu = ShiftMap.to(param.device)
        ShiftMap_gpu = torch.nn.Parameter(ShiftMap_gpu, requires_grad=True)
        ShiftMap_gpu2 = ShiftMap_gpu
        self.ShiftMap = ShiftMap_gpu
        self.ShiftMap2 = ShiftMap_gpu.clone().detach().requires_grad_(True)
        self.param = param

    def forward(self, Phantom):
        Mesh = Initialize_pathsmesh_interp_unwarp_2d(self.param.PhantomSize)
        Mesh = Mesh.astype('float32')
        Mesh = torch.tensor(Mesh).to(self.param.device)

        ShiftMap = self.ShiftMap

        ShiftMap_interp = PathsMeshingInterpolation_2d(ShiftMap, self.param)

        PhantomSize = self.param.PhantomSize
        Phantom_rot = Phantom.reshape(1, 1, PhantomSize[0], PhantomSize[1])

        ShiftMesh = ShiftMap_interp + Mesh

        warped_Phantom = F.grid_sample(Phantom_rot, ShiftMesh,
                                       mode='bicubic',
                                       padding_mode='zeros',
                                       align_corners=True)

        return warped_Phantom, ShiftMap_interp


def PathsMeshingInterpolation_2d(pathsmesh, param):
    pathsmesh = pathsmesh.permute(0, 3, 1, 2)
    pathsmesh_interp = F.interpolate(pathsmesh,
                                     size=tuple(param.PhantomSize),
                                     scale_factor=None,
                                     mode='bicubic')

    pathsmesh_interp = pathsmesh_interp.permute(0, 2, 3, 1)
    return pathsmesh_interp


def Initialize_ShiftMap_ept(Num_Ray):

    ShiftMap_x = np.zeros((Num_Ray[0], Num_Ray[1])).astype('float32') + 1e-8
    ShiftMap_y = np.zeros((Num_Ray[0], Num_Ray[1])).astype('float32') + 1e-8
    ShiftMap = np.stack([ShiftMap_x, ShiftMap_y], 2)

    return ShiftMap


def Initialize_pathsmesh_interp_unwarp_2d(PhantomSize):

    X_order = np.linspace(-1, 1, PhantomSize[1]).astype(np.float32)
    Y_order = np.linspace(-1, 1, PhantomSize[0]).astype(np.float32)

    [Phantom_X, Phantom_Y] = np.meshgrid(X_order, Y_order)

    mesh = np.stack([Phantom_X, Phantom_Y], 2)
    mesh = mesh[np.newaxis, :]
    return mesh


N = int(sys.argv[1])
M = int(sys.argv[2])
stop_epoch = int(sys.argv[3])
ii = int(sys.argv[4])
jj = int(sys.argv[5])
sidelobe = int(sys.argv[6])
lr = float(sys.argv[7])
Meta_name = sys.argv[8]
save_name = sys.argv[9]


Phantom_stack = tiff3Dread(Meta_name+'Meta_No'+'_'+str(ii)+'_'+str(jj)+'.tif')
Phantom_stack = np.transpose(Phantom_stack, [1, 2, 0])
Phantom_stack = Phantom_stack.astype('float32')
Phantom = torch.tensor(Phantom_stack[:, :, 0])

Phantom_fix = tiff3Dread(Meta_name+'Meta_No'+'_'+str(8)+'_'+str(8)+'.tif')
Phantom_fix = np.transpose(Phantom_fix, [1, 2, 0])
Phantom_fix = Phantom_fix.astype('float32')
Phantom_fix = torch.tensor(Phantom_fix[:, :, 0])

# Parameter
param = parameter()
device = torch.device("cuda:0")
param.device = device
param.PhantomSize = list(Phantom.shape)

# trans gpu
Phantom_gpu = Phantom.to(param.device)
Phantom_fix_gpu = Phantom_fix.reshape(1, 1, Phantom.shape[0], Phantom.shape[1]).to(param.device)

Model = Modle_A_G(param, N, M).to(param.device)
criterion = torch.nn.MSELoss(reduction='sum').to(device)# 交叉熵
optimizer = torch.optim.Adam([Model.ShiftMap], lr=lr)

##
ProjectName = save_name + 'Meta_'+str(ii)+'_'+str(jj)+'/'
SavePath_Map = ProjectName + str(N) + '_' + str(M) + '_' + str(lr) + '/Shift_map/'
SavePath_Img = ProjectName + str(N) + '_' + str(M) + '_' + str(lr) + '/Warped_Img/'
SavePath_Loss = ProjectName + str(N) + '_' + str(M) + '_' + str(lr) + '/'
if not os.path.exists(SavePath_Map):
    os.makedirs(SavePath_Map)
if not os.path.exists(SavePath_Img):
    os.makedirs(SavePath_Img)

# Training
Loss_list = np.zeros(stop_epoch)
lr_list = np.zeros(stop_epoch)
for epoch in range(0, stop_epoch):
    # 训练部分
    Model.train()
    Time_1 = time()

    warped_Phantom, ShiftMap_interp = Model(Phantom_gpu)

    # 计算损失函数
    loss = criterion(torch.squeeze(warped_Phantom[:, :, sidelobe:-1-sidelobe, sidelobe:-1-sidelobe]),
                     torch.squeeze(Phantom_fix_gpu[:, :, sidelobe:-1-sidelobe, sidelobe:-1-sidelobe]))

    optimizer.zero_grad()
    loss.backward()  # 背向传播
    optimizer.step()  # 优化器进行更新

    Loss_list[epoch] = loss.item()
    io.savemat(SavePath_Loss + str(lr) + '_Loss.mat', {'Loss': Loss_list})

    sys.stdout.write("[Train] [Epoch {}/{}] [loss:{:.8f}] time {:.3f}\n"
                     .format(epoch + 1, stop_epoch, loss.item(), time() - Time_1))
    sys.stdout.flush()

    if (epoch % 20 == 0) | (epoch == stop_epoch-1):
        tmp = np.squeeze(np.array(warped_Phantom.cpu().detach().numpy()))
        tmp = tmp / np.max(tmp)
        warped_Phantom_numpy = tmp * 65535
        warped_Phantom_numpy = warped_Phantom_numpy.astype(np.uint16)
        addr = SavePath_Img + 'warped_Phantom_' + str(epoch)
        io.savemat('addr', {'warped_Phantom': warped_Phantom_numpy})
        im = Image.fromarray(warped_Phantom_numpy)
        im.save(addr + '.tif')

        ShiftMap_gpu_numpy = np.squeeze(np.array(Model.ShiftMap.cpu().detach().numpy()))
        addr = SavePath_Map + 'ShiftMap_' + str(epoch) + '.mat'
        savemat(addr, {'ShiftMap': ShiftMap_gpu_numpy})
