from math import nan
from types import NoneType
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mplt

def scale_annotation(ax, xmax, ymin, offset, length, width, scale_text, text_offset, annotation_color, fontsize = 15):

    # Line 
    ax.arrow(xmax - length - offset, ymin + offset, length, 0.0, width = width, color = annotation_color,
        length_includes_head = True, head_width = 0.0)

    # x label
    text_x = xmax - 0.5 * length - offset
    text_y = ymin + text_offset

    ax.text(text_x, text_y, scale_text, color = annotation_color, fontsize = fontsize, horizontalalignment = "center", verticalalignment = "center")



def time_annotation(ax, xmin, ymax, offset, z_text, annotation_color):

    # x label
    text_x = xmin + offset
    text_y = ymax - offset

    ax.text(text_x, text_y, z_text, color = annotation_color, fontsize = 20, horizontalalignment = "left", verticalalignment = "top")



def text_annotation(ax, xmin, ymax, offset, z_text, annotation_color):

    # x label
    text_x = xmin - offset
    text_y = ymax - offset

    ax.text(text_x, text_y, z_text, color = annotation_color, fontsize = 20, horizontalalignment = "right", verticalalignment = "top")



"""
    propaganda_plot_columns(Nrows, Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name;
                            map_arr = nothing, par_arr = nothing,
                            contour_arr = nothing, contour_par_arr = nothing,
                            log_map = trues(Ncols),
                            colorbar_bottom = false,
                            colorbar_single = false,
                            smooth_col = falses(Ncols),
                            smooth_size = 0.0,
                            streamline_files = nothing,
                            streamlines = falses(Ncols),
                            contour_files = nothing,
                            contours = falses(Ncols),
                            contour_levels = nothing,
                            contour_color = "white",
                            smooth_contour_col = falses(Ncols),
                            alpha_contours = ones(Ncols),
                            cutoffs = nothing,
                            mask_bad = trues(Ncols),
                            bad_colors = ["k" for _ = 1:Ncols],
                            annotate_time = falses(Nrows * Ncols),
                            time_labels = nothing,
                            annotate_text = falses(Nrows * Ncols),
                            text_labels = nothing,
                            annotate_scale = trues(Nrows * Ncols),
                            annotation_color = "w",
                            scale_label = L"1 \\: h^{-1} c" * "Mpc",
                            scale_kpc = 1000.0,
                            r_circles = [0.0, 0.0],
                            shift_colorbar_labels_inward = trues(Ncols),
                            upscale = Ncols + 0.5,
                            read_mode = 1,
                            image_num = ones(Int64, Ncols),
                            transparent = false,
                            ticks_color = "k"
                            )

Creates an `image_grid` plot with `Ncols` and `Nrows` with colorbars on top.

"""

def propaganda_plot_columns(Nrows, Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name,\
                                map_arr = None, par_arr =  None,\
                                contour_arr =  None, contour_par_arr =  None,\
                                log_map = 'noentry',\
                                colorbar_bottom = False,\
                                colorbar_single = False,\
                                smooth_col = 'noentry',\
                                smooth_sizes = 0.0,\
                                annotate_smoothing = 'noentry',\
                                streamline_files =  None,\
                                streamlines = 'noentry',\
                                contour_files =  None,\
                                contours = 'noentry',\
                                contour_levels =  None,\
                                contour_color = "white",\
                                smooth_contour_col = 'noentry',\
                                alpha_contours = 'noentry',\
                                cutoffs =  None,\
                                mask_bad = 'noentry',\
                                bad_colors = 'noentry',\
                                annotate_time = 'noentry',\
                                time_labels =  None,\
                                annotate_text = 'noentry',\
                                text_labels =  None,\
                                annotate_scale = 'noentry',\
                                annotation_color = "w",\
                                scale_label = "1 \: h^{-1} c" * "Mpc",\
                                scale_kpc = 1000.0,\
                                r_circles = [0.0, 0.0],\
                                shift_colorbar_labels_inward = 'noentry',\
                                upscale = 'noentry',\
                                read_mode = 1,\
                                image_num = 'noentry',\
                                transparent = False,\
                                ticks_color = "k"\
                            ):

    if log_map == 'noentry':
        log_map = np.ones(Ncols, dtype=bool)
    if smooth_col == 'noentry':
         smooth_col=np.zeros(Ncols, dtype=bool)
    if annotate_smoothing == 'noentry':
        annotate_smoothing = np.zeros(Ncols, dtype=bool)
    if streamlines == 'noentry' :
        streamlines = np.zeros(Ncols, dtype=bool)
    if contours == 'noentry':
        contours = np.zeros(Ncols, dtype=bool)
    if smooth_contour_col == 'noentry':
        smooth_contour_col = np.zeros(Ncols, dtype=bool)
    if alpha_contours == 'noentry':
        alpha_contours = np.ones(Ncols)
    if mask_bad == 'noentry':
        mask_bad = np.ones(Ncols, dtype=bool)
    if  bad_colors == 'noentry':
        bad_colors = ["k"]*Ncols
    if annotate_time == 'noentry':
        annotate_time = np.zeros(Nrows * Ncols, dtype=bool)
    if annotate_text == 'noentry':
        annotate_text = np.zeros(Nrows * Ncols, dtype=bool)
    if annotate_scale == 'noentry':
        annotate_scale = np.ones(Nrows * Ncols, dtype=bool)
    if shift_colorbar_labels_inward == 'noentry':
        shift_colorbar_labels_inward = np.ones( Ncols, dtype=bool)
    if upscale == 'noentry':
        upscale = Ncols + 0.5
    if image_num == 'noentry':
        image_num = np.ones(Ncols)


    axes_grid1 = "mpl_toolkits.axes_grid1"

    axis_label_font_size = 20
    tick_label_size = 15

    fig = plt.figure(figsize = (6 * upscale * Nrows, 6.5 * upscale * Ncols))
    plt.plot_styling(color=ticks_color)

    # as in plt.subplot(111)
    grid = plt.axes_grid1.ImageGrid(fig, 111, \
    nrows_ncols = (Nrows, Ncols),\
    axes_pad = 0.0,\
    direction = "column",\
    cbar_location = "top",\
    cbar_mode = "edge",\
    cbar_size = "7%",\
    cbar_pad = 0.0)

    Nfile = 1
    Ncontour = 1
    time_label_num = 1
    
    for col in range(Ncols):

        for i in range(Nrows):

            print("Info: Column", col, "Plot", i)

            ax = grid[(col-1)*Nrows+i]

            if map_arr == None:

                image_name = files[Nfile]

                # read SPHtoGrid image
                if read_mode == 1:
                    map, par, snap_num, units = read_fits_image(image_name)

                    # read Smac1 binary image
                elif read_mode == 2:
                    map = read_smac1_binary_image(image_name)
                    smac1_info = read_smac1_binary_info(image_name)
                    smac1_center = [smac1_info.xcm, smac1_info.ycm, smac1_info.zcm]
                    par = mappingParameters(center = smac1_center,
                        x_size = smac1_info.boxsize_kpc,
                        y_size = smac1_info.boxsize_kpc,
                        z_size = smac1_info.boxsize_kpc,
                        Npixels = smac1_info.boxsize_pix)

                elif read_mode == 3:
                    map = read_smac2_image(image_name, image_num[Nfile])
                    par = read_smac2_info(image_name)
     
            else:
                map = map_arr[Nfile]
                par = par_arr[Nfile]
   

            if smooth_col[col]:
                pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
                smooth_pixel = float(smooth_sizes[col])/ pixelSideLength
                map = imfilter(map, reflect(Kernel.gaussian((smooth_pixel[1], smooth_pixel[2]), (par.Npixels[1] - 1, par.Npixels[1] - 1))))
  
            #!isnothing ????
            if cutoffs != None:
                for i in range(len(map)):
                    if map[i]<cutoffs[Nfile]:
                        map[i] = nan
   

            if mask_bad[col]:
                # get colormap object
                cmap = plt.get_cmap(im_cmap[col])
                # set invalid pixels to bad color
                cmap.set_bad(bad_colors[col])
            else:
                # get colormap object
                cmap = plt.get_cmap(im_cmap[col])
                # set invalid pixels to minimum
                cmap.set_bad(cmap(vmin_arr[col]))
   

            if log_map[col]:

                im = ax.imshow(map, norm = mplt.colors.LogNorm(),
                    vmin = vmin_arr[col], vmax = vmax_arr[col],
                    cmap = im_cmap[col],
                    origin = "lower"
                )
            else:
                im = ax.imshow(map,
                    vmin = vmin_arr[col], vmax = vmax_arr[col],
                    cmap = im_cmap[col],
                    origin = "lower"
                )
       

            if contours[col]:
                image_name = contour_files[Nfile]

                if contour_arr == None:
                    if read_mode == 1:
                        map, par, snap_num, units = read_fits_image(image_name)

                        # read Smac1 binary image
                    elif read_mode == 2:
                        map = read_smac1_binary_image(image_name)
                        smac1_info = read_smac1_binary_info(image_name)
                        smac1_center = [smac1_info.xcm, smac1_info.ycm, smac1_info.zcm] / 3.085678e21
                        par = mappingParameters(center = smac1_center,
                            x_size = smac1_info.boxsize_kpc,
                            y_size = smac1_info.boxsize_kpc,
                            z_size = smac1_info.boxsize_kpc,
                            Npixels = smac1_info.boxsize_pix)

                    elif read_mode == 3:

                        map = read_smac2_image(image_name, image_num[Nfile])
                        par = read_smac2_info(image_name)


                else:
                    map = contour_arr[Nfile]
                    par = contour_par_arr[Nfile]
               

                if smooth_contour_col[col]:
                    pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
                    smooth_pixel    = float(smooth_sizes[col]) / pixelSideLength
                    map = imfilter(map, reflect(Kernel.gaussian((smooth_pixel[1], smooth_pixel[2]), (par.Npixels[1] - 1, par.Npixels[1] - 1))))
                

                if contour_levels == None:
                    ax.contour(map, colors = contour_color, linewidth = 1.2, linestyle = "--", alpha = alpha_contours[col])
                else:
                    ax.contour(map, contour_levels, colors = contour_color, linewidth = 1.2, linestyle = "--", alpha = alpha_contours[col])
                

                Ncontour += 1
            

            if streamlines[col]:
                image_name = streamline_files[Nfile]
                vx, par, snap_num, units = read_fits_image(image_name)
                image_name = streamline_files[Nfile+1]
                vy, par, snap_num, units = read_fits_image(image_name)


                #should start from 1 or 0
                x_grid = np.arange(par.Npixels[1])
                y_grid = np.arange(par.Npixels[1])

                ax.streamplot(x_grid, y_grid,
                    vx, vy,
                    density = 2,
                    color = "white"
                )

                Ncontour += 1
            

            if col == 1:

                pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]


                # annotate_arrows(ax, 1.0, 1.0, 300.0/4.0, 500.0/4.0, 
                # 		x_arrow_text[i], y_arrow_text[i], 200.0/4.0)

                if annotate_scale[i]:
                    # annotate_scale(ax, par.Npixels[1], 1.0, 300.0/4.0, 1000.0/pixelSideLength, 
                    # 	L"1 \: h^{-1} c" * "Mpc", 500.0/4.0)

                    scale_annotation(ax, par.Npixels[1], 1.0, par.Npixels[1] / 14, \
                        scale_kpc / pixelSideLength, 0.1, \
                        scale_label, \
                        par.Npixels[1] / 8, \
                        annotation_color \
                    )
                

        

            if annotate_time[Nfile]:
                time_annotation(ax, 1.0, par.Npixels[1], 0.075 * par.Npixels[1],
                                time_labels[Nfile], annotation_color)
            

            if annotate_text[Nfile]:
                text_annotation(ax, par.Npixels[1], par.Npixels[1], 0.075 * par.Npixels[1],
                                text_labels[Nfile], annotation_color)
            

            pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]

            for r_circle in r_circles:
                ax.add_artist(plt.Circle((0.5*par.Npixels[1], 0.5*par.Npixels[2]), r_circle / pixelSideLength,\
                    color = "w", fill = False, ls = "--"))
           

            # draw smoothing beam
            if smooth_col[col] or smooth_contour_col[col] or annotate_smoothing[col]:
                pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
                smooth_pixel    = float(smooth_sizes[col]) / pixelSideLength
# anchor
# width
# height
                ax.add_artist(mplt.patches.Rectangle((0.1*par.Npixels[1] - 1.5*smooth_pixel[1],0.1*par.Npixels[1] - 1.5*smooth_pixel[2]),\
                                                            3*smooth_pixel[1], \
                                                            3*smooth_pixel[2], \
                                                            linewidth=1,edgecolor="gray",facecolor="gray"))

                ax.add_artist(mplt.patches.Ellipse((0.1*par.Npixels[1], 0.1*par.Npixels[2]), smooth_pixel[1], smooth_pixel[2], \
                    color = "w", fill = True, ls = "-"))
           

            ax.set_axis_off()

            if i == 1 or ((i == Nrows) and colorbar_bottom):

                if colorbar_single == True:
                    if (i == 1 or ((i == Nrows) and colorbar_bottom)):
                        Nfile += 1
                        continue
                    
                    #cax,kw = make_axes([grid[cax_i].cax for cax_i ∈ 1:Ncols])
                    cb = plt.colorbar(im, cax = grid[1], orientation = "horizontal")
                else:
                    cb = plt.colorbar(im, cax = grid[(col-1)*Nrows+1].cax, orientation = "horizontal")
               

                cb.set_label(cb_labels[col], fontsize = axis_label_font_size)
                cb.ax.tick_params( \
                    direction = "in",\
                    which = "major",\
                    labelsize = tick_label_size,\
                    size = 6, width = 1\
                )
                cb.ax.tick_params(\
                    direction = "in",\
                    which = "minor", \
                    labelsize = tick_label_size,\
                    size = 3, width = 1\
                )
                # cb_ticks_styling!(cb.ax, color=ticks_color)


                grid[(col-1)*Nrows+1].cax.xaxis.set_ticks_position("top")
                grid[(col-1)*Nrows+1].cax.xaxis.set_label_position("top")

                if shift_colorbar_labels_inward[col]:
                    plt.shift_colorbar_label(grid[(col-1)*Nrows+1].cax, "left")
                    plt.shift_colorbar_label(grid[(col-1)*Nrows+1].cax, "right")
                
            

            for spine in ax.spines.values():
                spine.set_edgecolor(ticks_color)
            

            # count up file
            Nfile += 1
        

    

    print("Info: saving", plot_name)
    plt.savefig(plot_name, bbox_inches = "tight", transparent=transparent)

    plt.close(fig)




"""
    propaganda_plot_columns(Nrows, Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name;
                            map_arr = nothing, par_arr = nothing,
                            contour_arr = nothing, contour_par_arr = nothing,
                            log_map = trues(Ncols),
                            colorbar_bottom = false,
                            colorbar_single = false,
                            smooth_col = falses(Ncols),
                            smooth_size = 0.0,
                            streamline_files = nothing,
                            streamlines = falses(Ncols),
                            contour_files = nothing,
                            contours = falses(Ncols),
                            contour_levels = nothing,
                            contour_color = "white",
                            smooth_contour_col = falses(Ncols),
                            alpha_contours = ones(Ncols),
                            cutoffs = nothing,
                            mask_bad = trues(Ncols),
                            bad_colors = ["k" for _ = 1:Ncols],
                            annotate_time = falses(Nrows * Ncols),
                            time_labels = nothing,
                            annotate_text = falses(Nrows * Ncols),
                            text_labels = nothing,
                            annotate_scale = trues(Nrows * Ncols),
                            annotation_color = "w",
                            scale_label = L"1 \\: h^{-1} c" * "Mpc",
                            scale_kpc = 1000.0,
                            r_circles = [0.0, 0.0],
                            shift_colorbar_labels_inward = trues(Ncols),
                            upscale = Ncols + 0.5,
                            read_mode = 1,
                            image_num = ones(Int64, Ncols),
                            transparent = false,
                            ticks_color = "k"
                            )

Creates an `image_grid` plot with `Ncols` and `Nrows` with colorbars on top.

"""

#poner el Int64
#scale_label = L"1 \: h^{-1} c" * "Mpc"
def propaganda_plot_rows(Nrows, Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name, \
                                map_arr = None, par_arr = None,\
                                contour_arr = None, contour_par_arr = None,\
                                log_map = 'noentry',\
                                colorbar_bottom = False,\
                                colorbar_single = False,\
                                smooth_row = 'noentry',\
                                smooth_sizes = 0.0,
                                annotate_smoothing = 'noentry', \
                                streamline_files = None, \
                                streamlines = 'noentry',\
                                contour_files = None,\
                                contours = 'noentry',\
                                contour_levels = None, \
                                contour_color = "white", \
                                smooth_contour_row = 'noentry', \
                                alpha_contours = 'noentry', \
                                cutoffs = None,  \
                                mask_bad = 'noentry', \
                                bad_colors = 'noentry', \
                                annotate_time = 'noentry',\
                                time_labels = None, \
                                annotate_text = 'noentry',\
                                text_labels = None, \
                                annotate_scale = 'noentry', \
                                annotation_color = "w", \
                                scale_label = "1 \: h^{-1} c" * "Mpc",\
                                scale_kpc = 1000.0,\
                                r_circles = [0.0, 0.0],\
                                shift_colorbar_labels_inward = 'noentry',\
                                upscale = 'noentry',\
                                read_mode = 1,\
                                image_num = 'noentry',\
                                transparent = False,\
                                ticks_color = "k"):



    if log_map == 'noentry':
        log_map = np.ones(Nrows, dtype=bool)
    if smooth_row == 'noentry':
        smooth_row = np.zeros(Nrows, dtype=bool)
    if annotate_smoothing ==  'noentry':
        annotate_smoothing = np.zeros(Nrows, dtype=bool)
    if streamlines ==  'noentry':
        streamlines = np.zeros(Nrows, dtype=bool)
    if contours ==  'noentry':
        contours = np.zeros(Nrows, dtype=bool)
    if smooth_contour_row ==  'noentry':
        smooth_contour_row = np.zeros(Nrows, dtype=bool)
    if alpha_contours ==  'noentry':
        alpha_contours = np.ones(Nrows)
    if mask_bad ==  'noentry':
        mask_bad = np.ones(Nrows, dtype=bool)
    if bad_colors ==  'noentry':
        bad_colors = ["k"]*Ncols
    if annotate_time ==  'noentry':
        annotate_time = np.zeros(Nrows* Ncols, dtype=bool)
    if annotate_text ==  'noentry':
        annotate_text = np.zeros(Nrows* Ncols, dtype=bool)
    if annotate_scale ==  'noentry':
        annotate_scale = np.ones(Ncols, dtype=bool)
    if shift_colorbar_labels_inward ==  'noentry':
        shift_colorbar_labels_inward = np.ones(Nrows, dtype=bool)
    if upscale ==  'noentry':
        upscale = Ncols + 0.5
    if image_num ==  'noentry':
        image_num = np.ones(Nrows)
    


    axes_grid1 = "mpl_toolkits.axes_grid1"

    axis_label_font_size = 20
    tick_label_size = 15

    fig = plt.figure(figsize = (6 * upscale * Nrows, 6.5 * upscale * Ncols))
    plt.plot_styling(color=ticks_color)

    grid = axes_grid1.ImageGrid(fig, 111,          # as in plt.subplot(111)
        nrows_ncols = (Nrows, Ncols),
        axes_pad = 0.0,
        direction = "column",
        cbar_location = "right",
        cbar_mode = "edge",
        cbar_size = "7%",
        cbar_pad = 0.0,
    )

    Nfile = 1
    Ncontour = 1
    time_label_num = 1
    
    for i in range(Nrows):
        for col in range(Ncols):

            print("Info: Column ", col, " Plot ",i)

            ax = grid[(col-1)*Nrows+i]

            if map_arr == None:

                image_name = files[Nfile]

                # read SPHtoGrid image
                if read_mode == 1:
                    map, par, snap_num, units = read_fits_image(image_name)

                    # read Smac1 binary image
                elif read_mode == 2:
                    map = read_smac1_binary_image(image_name)
                    smac1_info = read_smac1_binary_info(image_name)
                    smac1_center = [smac1_info.xcm, smac1_info.ycm, smac1_info.zcm]
                    par = mappingParameters(center = smac1_center,
                        x_size = smac1_info.boxsize_kpc,
                        y_size = smac1_info.boxsize_kpc,
                        z_size = smac1_info.boxsize_kpc,
                        Npixels = smac1_info.boxsize_pix)

                elif read_mode == 3:
                    map = read_smac2_image(image_name, image_num[Nfile])
                    par = read_smac2_info(image_name)
                
            else:
                map = map_arr[Nfile]
                par = par_arr[Nfile]
            

            if smooth_row[i]:
                pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
                smooth_pixel = float(smooth_sizes[i]) / pixelSideLength
                map = imfilter(map, reflect(Kernel.gaussian((smooth_pixel[1], smooth_pixel[2]), (par.Npixels[1] - 1, par.Npixels[1] - 1))))
            

            if cutoffs != None:
                for i in range(len(map)):
                    if map[i]<cutoffs[Nfile]:
                        map[i] = nan
            

            if mask_bad[i]:
                # get colormap object
                cmap = plt.get_cmap(im_cmap[i])
                # set invalid pixels to bad color
                cmap.set_bad(bad_colors[i])
            else:
                # get colormap object
                cmap = plt.get_cmap(im_cmap[i])
                # set invalid pixels to minimum
                cmap.set_bad(cmap(vmin_arr[i]))
            

            if log_map[i]:

                im = ax.imshow(map, norm = mplt.colors.LogNorm(),
                    vmin = vmin_arr[i], vmax = vmax_arr[i],
                    cmap = im_cmap[i],
                    origin = "lower"
                )
            else:
                im = ax.imshow(map,
                    vmin = vmin_arr[i], vmax = vmax_arr[i],
                    cmap = im_cmap[i],
                    origin = "lower"
                )
            

            if contours[i]:
                image_name = contour_files[Nfile]

                if contour_arr == None:
                    if read_mode == 1:
                        map, par, snap_num, units = read_fits_image(image_name)

                        # read Smac1 binary image
                    elif read_mode == 2:
                        map = read_smac1_binary_image(image_name)
                        smac1_info = read_smac1_binary_info(image_name)
                        smac1_center = [smac1_info.xcm, smac1_info.ycm, smac1_info.zcm] / 3.085678e21
                        par = mappingParameters(center = smac1_center,
                            x_size = smac1_info.boxsize_kpc,
                            y_size = smac1_info.boxsize_kpc,
                            z_size = smac1_info.boxsize_kpc,
                            Npixels = smac1_info.boxsize_pix)

                    elif read_mode == 3:

                        map = read_smac2_image(image_name, image_num[Nfile])
                        par = read_smac2_info(image_name)

                    
                else:
                    map = contour_arr[Nfile]
                    par = contour_par_arr[Nfile]
                

                if smooth_contour_row[i]:
                    pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
                    smooth_pixel    = smooth_sizes[i] / pixelSideLength
                    map = imfilter(map, reflect(Kernel.gaussian((smooth_pixel[1], smooth_pixel[2]), (par.Npixels[1] - 1, par.Npixels[1] - 1))))
                

                if contour_levels == None :
                    ax.contour(map, colors = contour_color, linewidth = 1.2, linestyle = "--", alpha = alpha_contours[i])
                else:
                    ax.contour(map, contour_levels, colors = contour_color, linewidth = 1.2, linestyle = "--", alpha = alpha_contours[i])
                

                Ncontour += 1
            

            if streamlines[i]:
                image_name = streamline_files[Nfile]
                vx, par, snap_num, units = read_fits_image(image_name)
                image_name = streamline_files[Nfile+1]
                vy, par, snap_num, units = read_fits_image(image_name)

                x_grid = np.arange(par.Npixels[1])
                y_grid = np.arange(par.Npixels[1])

                ax.streamplot(x_grid, y_grid, \
                    vx, vy,\
                    density = 2,\
                    color = "white" \
                )

                Ncontour += 1
            

            if col == 1 :

                pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]


                # annotate_arrows(ax, 1.0, 1.0, 300.0/4.0, 500.0/4.0, 
                # 		x_arrow_text[i], y_arrow_text[i], 200.0/4.0)

                if annotate_scale[i]:
                    # annotate_scale(ax, par.Npixels[1], 1.0, 300.0/4.0, 1000.0/pixelSideLength, 
                    # 	L"1 \: h^{-1} c" * "Mpc", 500.0/4.0)

                    scale_annotation(ax, par.Npixels[1], 1.0, par.Npixels[1] / 14, #scale_pixel_offset, 
                        scale_kpc / pixelSideLength, 0.1,
                        scale_label,
                        par.Npixels[1] / 8,
                        annotation_color
                    )
                

            

            if annotate_time[Nfile]:
                time_annotation(ax, 1.0, par.Npixels[1], 0.075 * par.Npixels[1],
                                time_labels[Nfile], annotation_color)
            

            if annotate_text[Nfile]:
                text_annotation(ax, par.Npixels[1], par.Npixels[1], 0.075 * par.Npixels[1],
                                text_labels[Nfile], annotation_color)
            

            pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]

            for r_circle in r_circles:
                ax.add_artist(plt.Circle((0.5*par.Npixels[1], 0.5*par.Npixels[2]), r_circle / pixelSideLength,
                    color = "w", fill = False, ls = "--"))
            

            # draw smoothing beam
            if smooth_row[i] or smooth_contour_row[i] or annotate_smoothing[i]:
                pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
                smooth_pixel    = smooth_sizes[i] / pixelSideLength

                ax.add_artist(mplt.patches.Rectangle((0.1*par.Npixels[1] - 1.5*smooth_pixel[1],0.1*par.Npixels[1] - 1.5*smooth_pixel[2]),           # anchor
                                                            3*smooth_pixel[1], # width
                                                            3*smooth_pixel[2], # height
                                                            linewidth=1,edgecolor="gray",facecolor="gray"))

                ax.add_artist(mplt.patches.Ellipse((0.1*par.Npixels[1], 0.1*par.Npixels[2]), smooth_pixel[1], smooth_pixel[2],
                    color = "w", fill = True, ls = "-"))
            

            ax.set_axis_off()

            if col == Ncols:

                if colorbar_single:
                    # if !(i == 1 || ((i == Nrows) && colorbar_bottom))
                    #     Nfile += 1
                    #     continue
                    # 
                    # #cax,kw = make_axes([grid[cax_i].cax for cax_i ∈ 1:Ncols])
                    # cb = colorbar(im, cax = grid[1])
                    pass
                else:
                    cb = plt.colorbar(im, cax = grid[i+Ncols-1].cax)
                

                cb.set_label(cb_labels[i], fontsize = axis_label_font_size)
                cb.ax.tick_params(
                    direction = "in",
                    which = "major",
                    labelsize = tick_label_size,
                    size = 6, width = 1
                )
                cb.ax.tick_params(
                    direction = "in",
                    which = "minor",
                    labelsize = tick_label_size,
                    size = 3, width = 1
                )
                # cb_ticks_styling!(cb.ax, color=ticks_color)


                if shift_colorbar_labels_inward[i]:
                    plt.shift_colorbar_label(grid[i+Ncols-1].cax, "bottom")
                    plt.shift_colorbar_label(grid[i+Ncols-1].cax, "top")
                
            

            for spine in ax.spines:
                spine.set_edgecolor(ticks_color)
            

            # count up file
            Nfile += 1
        

    

    print("Info: saving ", plot_name) 
    plt.savefig(plot_name, bbox_inches = "tight", transparent=transparent)

    plt.close(fig)








"""
    propaganda_plot_double_row(Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name;
                                map_arr = nothing, par_arr = nothing,
                                contour_arr = nothing, contour_par_arr = nothing,
                                log_map = trues(2Ncols),
                                smooth_file = falses(2Ncols),
                                smooth_sizes = 0.0,
                                annotate_smoothing = falses(2Ncols),
                                streamline_files = nothing,
                                streamlines = falses(2Ncols),
                                contour_files = nothing,
                                contours = falses(2Ncols),
                                contour_levels = nothing,
                                contour_color = "white",
                                smooth_contour_file = falses(2Ncols),
                                alpha_contours = ones(2Ncols),
                                cutoffs = nothing,
                                mask_bad = trues(2Ncols),
                                bad_colors = ["k" for _ = 1:2Ncols],
                                annotate_time = falses(2Ncols),
                                time_labels = nothing,
                                annotate_text = falses(2Ncols),
                                text_labels = nothing,
                                annotate_scale = trues(2Ncols),
                                scale_label = L"1 \\: h^{-1} c" * "Mpc",
                                scale_kpc = 1000.0,
                                r_circles = [0.0, 0.0],
                                shift_colorbar_labels_inward = trues(2Ncols),
                                upscale = Ncols,
                                aspect_ratio = 1.42,
                                read_mode = 1,
                                image_num = ones(Int64, 2Ncols),
                                transparent = false,
                                ticks_color = "k"
                                )

Creates an `image_grid` plot with `Ncols` and `Nrows` with colorbars on top.

"""


def propaganda_plot_double_row(Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name,
                                map_arr = None, par_arr = None,
                                contour_arr = None, contour_par_arr = None,
                                log_map = 'noentry',
                                smooth_file = 'noentry',
                                smooth_sizes = 0.0,
                                annotate_smoothing = 'noentry',
                                streamline_files =  None,
                                streamlines = 'noentry',
                                contour_files = None,
                                contours = 'noentry',
                                contour_levels = None,
                                contour_color = "white",
                                smooth_contour_file = 'noentry',
                                alpha_contours = 'noentry',
                                cutoffs = None,
                                mask_bad = 'noentry',
                                bad_colors = 'noentry',
                                annotate_time = 'noentry',
                                time_labels = None,
                                annotate_text = 'noentry',
                                text_labels = None,
                                annotate_scale = 'noentry',
                                scale_label = "1 \: h^{-1} c" * "Mpc",
                                scale_kpc = 1000.0,
                                r_circles = [0.0, 0.0],
                                shift_colorbar_labels_inward = 'noentry',
                                upscale ='noentry',
                                aspect_ratio = 1.42,
                                read_mode = 1,
                                image_num = 'noentry',
                                transparent = False,
                                ticks_color = "k"):

    if log_map == 'noentry':
        log_map = np.ones(2*Ncols, dtype=bool)
    if smooth_file == 'noentry':
        smooth_file = np.zeros(2*Ncols, dtype=bool)
    if annotate_smoothing == 'noentry':
        annotate_smoothing = np.zeros(2*Ncols, dtype=bool)
    if streamlines == 'noentry':
        streamlines = np.zeros(2*Ncols, dtype=bool)
    if contours == 'noentry':
        contours = np.zeros(2*Ncols, dtype=bool)
    if smooth_contour_file == 'noentry':
        smooth_contour_file = np.zeros(2*Ncols, dtype=bool)
    if alpha_contours == 'noentry':
        alpha_contours = np.ones(2*Ncols)
    if mask_bad == 'noentry':
        mask_bad = np.ones(2*Ncols, dtype=bool)
    if bad_colors == 'noentry':
        bad_colors = ["k"]* 2*Ncols -1
    if annotate_time == 'noentry':
        annotate_time = np.zeros(2*Ncols, dtype=bool)
    if annotate_text == 'noentry':
        annotate_text = np.zeros(2*Ncols, dtype=bool)
    if annotate_scale == 'noentry':
        annotate_scale = np.ones(2*Ncols, dtype=bool)
    if shift_colorbar_labels_inward == 'noentry':
        shift_colorbar_labels_inward = np.ones(2*Ncols, dtype=bool)
    if upscale == 'noentry':
        upscale = Ncols
    if image_num == 'noentry':
        image_num = np.ones(2*Ncols)

    fig = plt.get_figure( aspect_ratio , x_pixels = upscale*300)
    plt.plot_styling(upscale*300, axis_label_font_size=8, legend_font_size=5, color=ticks_color)

    gs = plt.GridSpec(4, Ncols, figure = fig, 
                      width_ratios = np.ones(Ncols), wspace = 0.0,
                      height_ratios = [0.05, 1, 1, 0.05], hspace=0.0
                      )

    Nfile = 1
    Ncontour = 1
    
    for row in range(2):
        for col in range(Ncols):

            print("Info: Column ", col, "Plot ",row)

            plt.subplot(get_gs(gs, row, (col-1)))
            ax = gca()

            if map_arr == None:

                image_name = files[Nfile]

                # read SPHtoGrid image
                if read_mode == 1:
                    map, par, snap_num, units = read_fits_image(image_name)

                    # read Smac1 binary image
                elif read_mode == 2:
                    map = read_smac1_binary_image(image_name)
                    smac1_info = read_smac1_binary_info(image_name)
                    smac1_center = [smac1_info.xcm, smac1_info.ycm, smac1_info.zcm] / 3.085678e21
                    par = mappingParameters(center = smac1_center,
                        x_size = smac1_info.boxsize_kpc,
                        y_size = smac1_info.boxsize_kpc,
                        z_size = smac1_info.boxsize_kpc,
                        Npixels = smac1_info.boxsize_pix)

                elif read_mode == 3:
                    map = read_smac2_image(image_name, image_num[Nfile])
                    par = read_smac2_info(image_name)
                
            else:
                map = map_arr[Nfile]
                par = par_arr[Nfile]
            

            if smooth_file[Nfile]:
                pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
                smooth_pixel    = smooth_sizes[Nfile] / pixelSideLength
                map = imfilter(map, reflect(Kernel.gaussian((smooth_pixel[1], smooth_pixel[2]), (par.Npixels[1] - 1, par.Npixels[1] - 1))))
            

            if cutoffs != None:

                for i in range(len(map)):
                    if map[i]<cutoffs[Nfile]:
                        map[i] = nan
            

            if mask_bad[Nfile] :
                # get colormap object
                cmap = plt.get_cmap(im_cmap[Nfile])
                # set invalid pixels to bad color
                cmap.set_bad(bad_colors[Nfile])
            

            if log_map[Nfile] :

                im = ax.imshow(map, norm = mplt.colors.LogNorm(),
                    vmin = vmin_arr[Nfile], vmax = vmax_arr[Nfile],
                    cmap = im_cmap[Nfile],
                    origin = "lower"
                )
            else:
                im = ax.imshow(map,
                    vmin = vmin_arr[Nfile], vmax = vmax_arr[Nfile],
                    cmap = im_cmap[Nfile],
                    origin = "lower"
                )
            

            if contours[Nfile]:
                image_name = contour_files[Nfile]

                if contour_arr == None:
                    if read_mode == 1:
                        map, par, snap_num, units = read_fits_image(image_name)

                        # read Smac1 binary image
                    elif read_mode == 2:
                        map = read_smac1_binary_image(image_name)
                        smac1_info = read_smac1_binary_info(image_name)
                        smac1_center = [smac1_info.xcm, smac1_info.ycm, smac1_info.zcm] / 3.085678e21
                        par = mappingParameters(center = smac1_center,
                            x_size = smac1_info.boxsize_kpc,
                            y_size = smac1_info.boxsize_kpc,
                            z_size = smac1_info.boxsize_kpc,
                            Npixels = smac1_info.boxsize_pix)

                    elif read_mode == 3:

                        map = read_smac2_image(image_name, image_num[Nfile])
                        par = read_smac2_info(image_name)

                    
                else:
                    map = contour_arr[Nfile]
                    par = contour_par_arr[Nfile]
                

                if smooth_contour_file[Nfile]:
                    #println("smooth size = $smooth_size kpc")
                    pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
                    smooth_pixel    = smooth_sizes[Nfile] / pixelSideLength
                    #println("smooth size = $smooth_size pix")
                    map = imfilter(map, reflect(Kernel.gaussian((smooth_pixel[1], smooth_pixel[2]), (par.Npixels[1] - 1, par.Npixels[1] - 1))))
                

                if contour_levels == None:
                    ax.contour(map, colors = contour_color, linewidth = 1.2, linestyle = "--", alpha = alpha_contours[col])
                else:
                    ax.contour(map, contour_levels, colors = contour_color, linewidth = 1.2, linestyle = "--", alpha = alpha_contours[col])
                

                Ncontour += 1
            

            if streamlines[Nfile]:
                image_name = streamline_files[Nfile]
                vx, par, snap_num, units = read_fits_image(image_name)
                image_name = streamline_files[Nfile+1]
                vy, par, snap_num, units = read_fits_image(image_name)

                x_grid = np.arange(par.Npixels[1])
                y_grid = np.arange(par.Npixels[1])

                ax.streamplot(x_grid, y_grid,\
                    vx, vy,\
                    density = 2,\
                    color = "white"\
                )

                Ncontour += 1
            

            if col == 1 :

                pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]


                # annotate_arrows(ax, 1.0, 1.0, 300.0/4.0, 500.0/4.0, 
                # 		x_arrow_text[i], y_arrow_text[i], 200.0/4.0)

                if annotate_scale[row]:
                    # annotate_scale(ax, par.Npixels[1], 1.0, 300.0/4.0, 1000.0/pixelSideLength, 
                    # 	L"1 \: h^{-1} c" * "Mpc", 500.0/4.0)

                    scale_annotation(ax, par.Npixels[1], 1.0, par.Npixels[1] / 14, #scale_pixel_offset, 
                        scale_kpc / pixelSideLength, 0.1,
                        scale_label,
                        par.Npixels[1] / 8,
                        annotation_color
                    )
                

            

            if annotate_time[Nfile]:
                time_annotation(ax, 1.0, par.Npixels[1], 0.075 * par.Npixels[1],
                                time_labels[Nfile],
                                annotation_color)
            

            if annotate_text[Nfile]:
                text_annotation(ax, par.Npixels[1], par.Npixels[1], 0.075 * par.Npixels[1],
                                text_labels[Nfile],
                                annotation_color)
            

            pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]

            for r_circle in r_circles:
                ax.add_artist(plt.Circle((0.5*par.Npixels[1], 0.5*par.Npixels[2]), r_circle / pixelSideLength,\
                    color = "w", fill = False, ls = "--"))
            

            # draw smoothing beam
            if smooth_file[Nfile] or smooth_contour_file[Nfile] or annotate_smoothing[Nfile]:
                pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
                smooth_pixel    = smooth_sizes[Nfile] / pixelSideLength

                ax.add_artist(mplt.patches.Rectangle((0.1*par.Npixels[1] - 1.5*smooth_pixel[1],0.1*par.Npixels[1] - 1.5*smooth_pixel[2]),           # anchor
                                                            3*smooth_pixel[1], # width
                                                            3*smooth_pixel[2], # height
                                                            linewidth=1,edgecolor="gray",facecolor="gray"))

                ax.add_artist(mplt.patches.Ellipse((0.1*par.Npixels[1], 0.1*par.Npixels[2]), smooth_pixel[1], smooth_pixel[2],
                    color = "w", fill = True, ls = "-"))
            

            ax.set_axis_off()

            for spine in ax.spines:
                spine.set_edgecolor(ticks_color)
            

            if row == 1:
                loc = "top"
                cax_row = 0
            else:
                loc = "bottom"
                cax_row = 3
            

            plt.subplot(get_gs(gs, cax_row, (col-1)))
            cax = gca()

            if log_map[Nfile]:
                sm = plt.cm.ScalarMappable(cmap = im_cmap[Nfile], norm = mplt.colors.LogNorm(vmin = vmin_arr[Nfile], vmax = vmax_arr[Nfile]))
            else:
                sm = plt.cm.ScalarMappable(cmap = im_cmap[Nfile], plt.Normalize(vmin = vmin_arr[Nfile], vmax = vmax_arr[Nfile]) )
            

            sm.set_array([])
            cb = plt.colorbar(sm, cax=cax, fraction = 0.046, orientation = "horizontal")

            cb.set_label(cb_labels[Nfile])
            cb.ax.tick_params(
                direction = "in",
                which = "major",
                size = 6, width = 1
            )
            cb.ax.tick_params(
                direction = "in",
                which = "minor",
                size = 3, width = 1
            )

            #cb_ticks_styling!(cb.ax, color=ticks_color)

            cb.ax.xaxis.set_ticks_position(loc)
            cb.ax.xaxis.set_label_position(loc)

            if shift_colorbar_labels_inward[Nfile]:
                shift_colorbar_label(cb.ax, "left")
                shift_colorbar_label(cb.ax, "right")
            
            

            # count up file
            Nfile += 1
        

    #gs.tight_layout(fig, pad = 0.0, h_pad = -10.0, w_pad= 0.0)

    print( "Info saving ", plot_name)
    plt.savefig(plot_name, bbox_inches = "tight", transparent=transparent)

    plt.close(fig)

