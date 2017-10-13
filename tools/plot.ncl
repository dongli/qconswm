load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
    print("[Error]: Please set file_prefix in command arguments!")
    exit
end if

if (isdefined("time_range")) then

else
    file_paths = systemfunc("ls " + file_prefix + "*.nc")
end if

wks = gsn_open_wks("pdf", file_prefix)

gsn_define_colormap(wks, "BlAqGrYeOrRe")

res = True
res@gsnDraw = False
res@gsnFrame = False

res_gd = res
res_gd@cnLineLabelsOn = False
res_gd@cnLevelSelectionMode = "ManualLevels"
res_gd@cnMinLevelValF = 7.5e4
res_gd@cnMaxLevelValF = 9.5e4
res_gd@cnLevelSpacingF = 0.1e4
res_gd@gsnPolar = "NH"

res_u = res

res_v = res

res_panel = True

do k = 0, dimsizes(file_paths) - 1
    print("[Notice]: Plot " + file_paths(k))
    f = addfile(file_paths(k), "r")
    time = f->time
    time@units = "hours since 0000-01-01 00:00:00"
    res_gd@gsnRightString = cd_string(f->time, "%J %H:%M")
    plot0 = gsn_csm_contour_map(wks, f->gd(0,:,:), res_gd)
    draw(plot0)
    frame(wks)
    ;plot1 = gsn_csm_contour_map(wks, f->u(0,:,:), res_u)
    ;plot2 = gsn_csm_contour_map(wks, f->v(0,:,:), res_v)
    ;gsn_panel(wks, (/plot0, plot1, plot2/), (/3, 1/), res_panel)
end do