load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

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
res@cnLinesOn = False
res@cnFillOn = True
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 7.5e4
res@cnMaxLevelValF = 9.5e4
res@cnLevelSpacingF = 0.1e4

do k = 0, dimsizes(file_paths) - 1
    print("[Notice]: Plot " + file_paths(k))
    f = addfile(file_paths(k), "r")
    plot = gsn_csm_contour_map(wks, f->gd(0,:,:), res)
end do