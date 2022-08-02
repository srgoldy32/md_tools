from Dump import Dump

class DumpSeries:
    def __init__(self, file_0, file_end, dump_step, ):
        self.file_0 = file_0
        self.file_end = file_end
        self.dump_step = dump_step
        
        self.d_0 = Dump(file_0)
        self.d_end = Dump(file_end)
        
        self.planes = self.get_planes()
        self.diffs = self.get_z_disp()
        self.RMS = self.get_RMS()
        
    def __str__(self):
        return "<Dump Series from: {} to {}>".format(self.file_0,self.file_end)
    
    def closest(self,lis, K):
        return lis[min(range(len(lis)), key = lambda i: abs(lis[i]-K))]

    def averaged_RMS(self,w):
        running_sum = 0
        N = len(w)
        for w_i in w:
            running_sum += w_i**2/N

        return running_sum**.5
    
    def get_planes(self):
        z_dict = {}
        for z in self.d_0.zs:
            if z in z_dict:
                z_dict[z] += 1
            else:
                z_dict[z] = 1
        
        self.planes = list(z_dict.keys())
        return self.planes
    
    def get_z_disp(self):
        self.z_disp = []
        for z in self.d_end.zs:
            home_plane = self.closest(self.planes,z)
            diff = abs(home_plane - z)
            self.z_disp.append(diff)
        
        return self.z_disp
    
    def get_RMS(self):
        return self.averaged_RMS(self.z_disp)
        