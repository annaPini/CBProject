from parameters import *
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button

# //////////////////////////////////////////////////////////////////////////////
class RMSDPlotter:
    init_ref_frame = 0
    init_rmsd_treshold = 2

    def __init__(self, rmsd_mat_path):
        ########## Generic
        print(f">>> Loading RMSD data...")
        self.rmsd_mat = np.load(rmsd_mat_path)

        ########## RMSD 1D
        self.x = np.arange(self.rmsd_mat.shape[0])


    def rmsd_2d(self):
        plt.colorbar(plt.imshow(self.rmsd_mat, vmax = 5))


    def rmsd_1d_sliders(self, WAD_button = False):
        plt.subplots()
        plt.subplots_adjust(bottom = 0.25, top = 0.9)

        self.rmsd_1d_color(self.init_ref_frame, self.init_rmsd_treshold)

        self.line_rmsd = plt.scatter(self.x, self.rmsd_arr, c = self.colors, marker = '.')
        self.line_treshold = plt.plot(self.x, self.y)

        plt.xlabel("Frame")
        plt.ylabel("RMSD")

        self.slid_ref_frame = Slider(ax = plt.axes((.15, .10, .5, .03)), label = "ref_frame", valmin = 0, valmax = self.rmsd_mat.shape[0] - 1,  valinit = self.init_ref_frame, valstep = 1,  color = "orange" )
        self.slid_rmsd_treshold = Slider(ax = plt.axes((.15, .05, .5, .03)), label = "rmsd_treshold", valmin = 0, valmax = np.max(self.rmsd_mat), valinit = self.init_rmsd_treshold, valstep = .05,  color = "blue" )

        self.slid_ref_frame.on_changed(self.rmsd_1d_update_plot)
        self.slid_rmsd_treshold.on_changed(self.rmsd_1d_update_plot)

        if WAD_button:
            self.button_save = Button(ax = plt.axes((.8, .05, .2, .1)), label = "Select WA frames")
            self.button_save.on_clicked(self.rmsd_1d_WAD_select_frames)


    def rmsd_1d_color(self, ref_frame, rmsd_treshold):
        frames = self.rmsd_mat.shape[0]
        self.rmsd_arr = self.rmsd_mat[ref_frame]

        self.colors = np.zeros((frames,3))
        self.colors[:,2] = 1
        self.colors[self.rmsd_arr <= rmsd_treshold] = GREEN
        self.colors[ref_frame] = RED

        self.y = np.zeros(frames) + rmsd_treshold


    def rmsd_1d_update_plot(self, val):
        self.rmsd_1d_color(self.slid_ref_frame.val, self.slid_rmsd_treshold.val)

        data = np.append(self.x, self.rmsd_arr).reshape((2, self.x.size)).T
        self.line_rmsd.set_offsets(data)
        self.line_rmsd.set_color(self.colors)

        self.line_treshold[0].set_data(self.x, self.y)

    def rmsd_1d_WAD_select_frames(self, val):
        ref_frame = self.slid_ref_frame.val
        rmsd_treshold = self.slid_rmsd_treshold.val

        self.rmsd_arr = self.rmsd_mat[ref_frame]
        frames = [int(f) for f in self.x[self.rmsd_arr <= rmsd_treshold]]

        info = Info(DIR_DA_WAD / f"{CURRENT_RUN}-{ref_frame}-info.json")

        info.update(
            rsmd_treshold = rmsd_treshold,
            ref_frame = ref_frame,
            ref_frame_index = frames.index(ref_frame),
            frames = frames,
        )

# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    for run in RUNS[:1]:
        PATH_RMSD = DIR_DA_GENERAL / f"{run}-rmsd.npy"

        plotter = RMSDPlotter(PATH_RMSD)
        print(f">>> Plotting RMSD data...")

        ###################################################################
        plotter.rmsd_1d_sliders(WAD_button = True)
        # plotter.rmsd_2d()


        plt.show()

# //////////////////////////////////////////////////////////////////////////////
