/* DSC Viewer V1.0.0: local NIfTI parsing and reproducible signal processing. */
(function(root){
'use strict';
const VERSION='1.0.0',MAX_VALUES=80000000;
function fail(s){throw new Error(s);}
function header(buffer){
 if(buffer.byteLength<352)fail('文件过短，无法读取 NIfTI 头。');
 const d=new DataView(buffer);let le=d.getInt32(0,true)===348;
 if(!le&&d.getInt32(0,false)!==348)fail('V1 支持 NIfTI-1 单文件 .nii/.nii.gz；该文件不是可识别的 NIfTI-1。DICOM、NIfTI-2、图片和 .hdr/.img 请先转换。');
 const i=o=>d.getInt16(o,le),f=o=>d.getFloat32(o,le);
 const magic=String.fromCharCode(...new Uint8Array(buffer,344,4));if(magic!=='n+1\0')fail('需要单文件 NIfTI（n+1），不支持独立 .hdr/.img。');
 const nd=i(40);if(nd<3||nd>7)fail('NIfTI 维数无效。');const dims=Array.from({length:7},(_,k)=>i(42+k*2));
 if(dims.slice(0,nd).some(x=>!Number.isInteger(x)||x<1))fail('NIfTI 尺寸包含非正整数。');
 if(dims.slice(4,nd).some(v=>v>1))fail('V1 不支持多回波/向量等五维数据，请先选择单个四维序列。');
 const [X,Y,Z]=dims,T=nd>=4?dims[3]:1,nvox=X*Y*Z,total=nvox*T;
 if(!Number.isSafeInteger(total)||total>MAX_VALUES)fail(`数据超过 V1 的 ${MAX_VALUES/1e6} 百万数值内存上限，请裁剪空间范围或拆分序列。`);
 const types={2:['getUint8',1],4:['getInt16',2],8:['getInt32',4],16:['getFloat32',4],64:['getFloat64',8],256:['getInt8',1],512:['getUint16',2],768:['getUint32',4]};
 const type=i(70),spec=types[type];if(!spec)fail('不支持此 NIfTI 数据类型：'+type+'。请转换为标量 int16 或 float32。');
 if(i(72)!==spec[1]*8)fail('数据类型与 bitpix 不一致。');
 const offset=f(108);if(!Number.isInteger(offset)||offset<352||offset+total*spec[1]>buffer.byteLength)fail('文件数据长度或 vox_offset 无效，文件可能不完整。');
 const pix=Array.from({length:8},(_,k)=>f(76+k*4));if(pix.slice(1,4).some(v=>!Number.isFinite(v)||v<=0))fail('空间间距无效，无法可靠显示物理距离。');
 const units=d.getUint8(123),spaceUnit=units&7,timeUnit=units&56;
 const sm={1:1000,2:1,3:.001}[spaceUnit]||1,tm={8:1,16:.001,24:.000001}[timeUnit];
 const warnings=[];if(![1,2,3].includes(spaceUnit))warnings.push('空间单位未声明：暂按毫米，需核实。');
 let dt=tm&&pix[4]>0&&Number.isFinite(pix[4])?pix[4]*tm:null;if(!dt)warnings.push('时间单位或 TR 未可靠声明：请确认界面的帧间隔。');
 let affine,affineSource;
 if(i(254)>0){affine=Array.from({length:3},(_,r)=>Array.from({length:4},(_,c)=>f(280+r*16+c*4)*sm));affineSource='sform';}
 else if(i(252)>0){let b=f(256),c=f(260),e=f(264),norm=b*b+c*c+e*e,a;
  if(norm>1){const q=Math.sqrt(norm);b/=q;c/=q;e/=q;a=0;}else a=Math.sqrt(1-norm);
  const sx=pix[1]*sm,sy=pix[2]*sm,sz=pix[3]*sm*(pix[0]<0?-1:1);
  affine=[[(a*a+b*b-c*c-e*e)*sx,2*(b*c-a*e)*sy,2*(b*e+a*c)*sz,f(268)*sm],[2*(b*c+a*e)*sx,(a*a+c*c-b*b-e*e)*sy,2*(c*e-a*b)*sz,f(272)*sm],[2*(b*e-a*c)*sx,2*(c*e+a*b)*sy,(a*a+e*e-c*c-b*b)*sz,f(276)*sm]];affineSource='qform';}
 else {affine=[[pix[1]*sm,0,0,0],[0,pix[2]*sm,0,0],[0,0,pix[3]*sm,0]];affineSource='pixdim';warnings.push('无有效 sform/qform：可显示采集索引，不能确认解剖方向。');}
 if(affine.flat().some(x=>!Number.isFinite(x)))fail('空间变换包含非法数值。');
 const spacing=[0,1,2].map(c=>Math.hypot(...affine.map(row=>row[c])));if(spacing.some(x=>x<=0))fail('空间变换退化。');
 let slope=f(112),inter=f(116);if(slope===0||!Number.isFinite(slope)){slope=1;inter=0;}if(!Number.isFinite(inter))fail('无效 scl_inter。');
 return {X,Y,Z,T,nvox,total,le,type,spec,offset,slope,inter,spacing,affine,affineSource,dt,units,warnings,version:VERSION};
}
function parse(buffer){const m=header(buffer),d=new DataView(buffer),values=new Float32Array(m.total);for(let i=0;i<m.total;i++)values[i]=d[m.spec[0]](m.offset+i*m.spec[1],m.le)*m.slope+m.inter;return {meta:m,values};}
function compatible(a,b){return a.X===b.X&&a.Y===b.Y&&a.Z===b.Z&&a.affine.every((r,i)=>r.every((v,j)=>Math.abs(v-b.affine[i][j])<1e-3));}
function quantile(a,q){const s=Array.from(a).filter(Number.isFinite).sort((a,b)=>a-b);if(!s.length)return 0;return s[Math.min(s.length-1,Math.floor((s.length-1)*q))];}
function kernel(s){if(!(s>0))return {r:0,w:[1]};const r=Math.ceil(3*s);if(r>64)fail('平滑尺度相对采样间隔过大，请减小平滑强度。');const w=[];let sum=0;for(let k=-r;k<=r;k++){const a=Math.exp(-k*k/(2*s*s));w.push(a);sum+=a;}return {r,w:w.map(v=>v/sum)};}
function temporal(raw,T,sigma){if(!(sigma>0))return raw.slice();const {r,w}=kernel(sigma),o=new Float32Array(raw.length);for(let p=0;p<raw.length/T;p++)for(let t=0;t<T;t++){let v=0;for(let k=-r;k<=r;k++)v+=w[k+r]*raw[p*T+Math.max(0,Math.min(T-1,t+k))];o[p*T+t]=v;}return o;}
function spatial(input,mask,m,sx,sy){
 const {X,Y,Z,T}=m;let weights=Float32Array.from(mask),v=input.map((a,i)=>mask[Math.floor(i/T)]?a:0);
 for(const [axis,s] of [[0,sx],[1,sy]]){if(!(s>0))continue;const {r,w}=kernel(s),next=new Float32Array(v.length),nw=new Float32Array(mask.length);
  for(let z=0;z<Z;z++)for(let y=0;y<Y;y++)for(let x=0;x<X;x++){const p=(z*Y+y)*X+x;for(let k=-r;k<=r;k++){const xx=x+(axis===0?k:0),yy=y+(axis===1?k:0);if(xx<0||xx>=X||yy<0||yy>=Y)continue;const q=(z*Y+yy)*X+xx,a=w[k+r];nw[p]+=a*weights[q];for(let t=0;t<T;t++)next[p*T+t]+=a*v[q*T+t];}}
  v=next;weights=nw;
 }
 for(let p=0;p<mask.length;p++)for(let t=0;t<T;t++)v[p*T+t]=mask[p]&&weights[p]>0?v[p*T+t]/weights[p]:0;
 return v;
}
function prepare(scan,cfg,externalMask,progress=()=>{}){
 const m=scan.meta,{X,Y,Z,T,nvox}=m;if(T<4)fail('需要至少4个时间点的四维动态序列，当前文件不能演示 DSC 动态。');
 const baseline=Number(cfg.baseline),dt=Number(cfg.dt),ts=Number(cfg.temporalSeconds),ss=Number(cfg.spatialMM),frac=Number(cfg.thresholdPercent)/100;
 if(!Number.isInteger(baseline)||baseline<1||baseline>T-2)fail('基线帧数必须为1至总帧数−2的整数。');
 if(!(dt>0)||!Number.isFinite(dt))fail('帧间隔必须为正数（秒）。');
 if(![ts,ss,frac].every(Number.isFinite)||ts<0||ss<0||frac<0||frac>1)fail('平滑或前景参数无效。');
 kernel(ts/dt);kernel(ss/m.spacing[0]);kernel(ss/m.spacing[1]);
 const bg=new Float32Array(nvox),finite=new Uint8Array(nvox),samples=[];
 for(let p=0;p<nvox;p++){let sum=0,ok=true;for(let t=0;t<T;t++){const v=scan.values[t*nvox+p];if(!Number.isFinite(v)||v<=0)ok=false;if(t<baseline)sum+=v;}bg[p]=sum/baseline;finite[p]=ok?1:0;if(Number.isFinite(bg[p])&&bg[p]>0&&p%Math.max(1,Math.floor(nvox/100000))===0)samples.push(bg[p]);}
 const p99=quantile(samples,.99);if(!(p99>0))fail('没有可用的正值信号。请确认输入是原始 DSC 幅度图像，而不是浓度、Tmax或其他参数图。');
 const threshold=p99*frac,mask=new Uint8Array(nvox),raw=new Float32Array(nvox*T);let count=0,excluded=0;
 for(let p=0;p<nvox;p++){let support=externalMask?externalMask[p]>0:bg[p]>threshold;if(support&&!finite[p])excluded++;if(!support||!finite[p]||bg[p]<=0)continue;mask[p]=1;count++;for(let t=0;t<T;t++)raw[p*T+t]=-Math.log(scan.values[t*nvox+p]/bg[p]);}
 if(!count)fail('当前阈值或掩膜下没有有效前景体素，请降低阈值或核对掩膜。');
 progress('时间平滑');let filtered=temporal(raw,T,ts/dt);if(ss>0){progress('切片内空间平滑');filtered=spatial(filtered,mask,m,ss/m.spacing[0],ss/m.spacing[1]);}else for(let p=0;p<nvox;p++)if(!mask[p])filtered.fill(0,p*T,(p+1)*T);
 const scale=quantile(samples,.995),gray=new Uint8Array(nvox);for(let p=0;p<nvox;p++)gray[p]=Math.round(Math.max(0,Math.min(1,bg[p]/Math.max(scale,1e-8)))*255);
 return {meta:{...m,dt},raw,filtered,mask,gray,audit:{version:VERSION,shape:[X,Y,Z,T],baselineFrames:baseline,dtSeconds:dt,timeSigmaSeconds:ts,timeSigmaFrames:ts/dt,spatialSigmaMM:ss,spatialSigmaVoxels:[ss/m.spacing[0],ss/m.spacing[1],0],foreground:externalMask?'external_mask':'baseline_signal_threshold_NOT_brain_extraction',thresholdPercent:cfg.thresholdPercent,thresholdAbsolute:threshold,validVoxels:count,excludedNonfiniteOrNonpositiveVoxels:excluded,formula:'per voxel -ln(S/mean(first baselineFrames S)); signed; no TE scaling; no extra centering',temporalBoundary:'nearest; normalized Gaussian radius ceil(3*sigma)',spatialMethod:'in-plane mask-normalized Gaussian; no between-slice smoothing or interpolation',outsideMask:'zero sentinel; not a measured zero',maskExcludesAnyNonpositiveOrNonfiniteTimeSample:true,grayWindow:[0,scale],rawSourceUnmodified:true}};
}
function jet(v,lo=-.1,hi=.65){const u=Math.max(0,Math.min(1,(v-lo)/(hi-lo)));const lerp=(xs,ys)=>{let i=1;while(i<xs.length-1&&u>xs[i])i++;return ys[i-1]+(ys[i]-ys[i-1])*(u-xs[i-1])/(xs[i]-xs[i-1]);};return [lerp([0,.35,.66,.89,1],[0,0,1,1,.5]),lerp([0,.125,.375,.64,.91,1],[0,0,1,1,0,0]),lerp([0,.11,.34,.65,1],[.5,1,1,0,0])].map(x=>Math.round(x*255));}
const api={VERSION,header,parse,prepare,compatible,kernel,temporal,spatial,jet,quantile};root.DSCV1=api;if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(globalThis);
