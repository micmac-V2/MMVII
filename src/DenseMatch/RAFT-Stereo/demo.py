import sys
#sys.path.append('core')
import time
import argparse
import glob
import numpy as np
import torch
from tqdm import tqdm
from pathlib import Path
from core.raft_stereo import RAFTStereo
from core.utils.utils import InputPadder
import tifffile as tf




def load_image(imfile, which_device):
    img = tf.imread(imfile).astype(np.uint8)
    if img.ndim==2:
        img = np.repeat(img[:, :, np.newaxis], 3, axis=2)
    img = torch.from_numpy(img).permute(2, 0, 1)
    return img[None].to(which_device)

def demo(args):
    DEVICE = torch.device('cuda:{}'.format(args.ongpu))
    model = torch.nn.DataParallel(RAFTStereo(args), device_ids=[args.ongpu])
    model.load_state_dict(torch.load(args.restore_ckpt,map_location=torch.device('cpu')))

    model = model.module
    model.to(DEVICE)
    model.eval()

    output_file = Path(args.output_file)

    with torch.no_grad():
        left_images = sorted(glob.glob(args.left_imgs, recursive=True))
        right_images = sorted(glob.glob(args.right_imgs, recursive=True))
        for (imfile1, imfile2) in tqdm(list(zip(left_images, right_images))):
            image1 = load_image(imfile1,DEVICE)
            image2 = load_image(imfile2,DEVICE)
            _, flow_up = model(image1, image2, iters=args.valid_iters, test_mode=True)
            tf.imwrite(output_file, flow_up.cpu().numpy().squeeze())



# duplicate demo to use torch profiler and checks bottlenecks
def demo_profile(args, n_warmup=3):
    from torch.profiler import profile, ProfilerActivity

    torch.backends.cudnn.benchmark = True
    DEVICE = torch.device(f"cuda:{args.ongpu}")
    model = torch.nn.DataParallel(RAFTStereo(args), device_ids=[args.ongpu])
    model.load_state_dict(torch.load(args.restore_ckpt, map_location="cpu"))
    model = model.module.to(DEVICE).eval()

    left = sorted(glob.glob(args.left_imgs, recursive=True))
    right = sorted(glob.glob(args.right_imgs, recursive=True))
    image1 = load_image(left[0], DEVICE)
    image2 = load_image(right[0], DEVICE)

    with torch.inference_mode():
        for _ in range(n_warmup):
            _, flow_up = model(image1, image2, iters=args.valid_iters, test_mode=True)
            del flow_up
        torch.cuda.synchronize()

        # 1) wall-clock timing, no profiler
        times = []
        for _ in range(10):
            torch.cuda.synchronize(); t0 = time.perf_counter()
            _, flow_up = model(image1, image2, iters=args.valid_iters, test_mode=True)
            torch.cuda.synchronize(); times.append((time.perf_counter() - t0) * 1000)
            del flow_up
        times.sort()
        print(f"median {times[len(times)//2]:.1f} ms, min {times[0]:.1f}, max {times[-1]:.1f}")
        print(f"peak mem {torch.cuda.max_memory_allocated()/2**30:.2f} GiB")

        # 2) one profiled pass
        with profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA]) as prof:
            _, flow_up = model(image1, image2, iters=args.valid_iters, test_mode=True)
            torch.cuda.synchronize()
        del flow_up

    print(prof.key_averages().table(sort_by="self_cuda_time_total", row_limit=25))
    prof.export_chrome_trace("raft_trace.json")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--restore_ckpt', help="restore checkpoint", required=True)
    parser.add_argument('--save_numpy', action='store_true', help='save output as numpy arrays')
    parser.add_argument('-l', '--left_imgs', help="path to all first (left) frames", default="datasets/Middlebury/MiddEval3/testH/*/im0.png")
    parser.add_argument('-r', '--right_imgs', help="path to all second (right) frames", default="datasets/Middlebury/MiddEval3/testH/*/im1.png")
    parser.add_argument('-o','--output_file', help="path to output image file", default="Px_.tif")
    parser.add_argument('--mixed_precision', action='store_true', help='use mixed precision')
    parser.add_argument('--valid_iters', type=int, default=32, help='number of flow-field updates during forward pass')

    # Architecture choices
    parser.add_argument('--hidden_dims', nargs='+', type=int, default=[128]*3, help="hidden state and context dimensions")
    parser.add_argument('--corr_implementation', choices=["reg", "alt", "reg_cuda", "alt_cuda"], default="reg", help="correlation volume implementation")
    parser.add_argument('--shared_backbone', action='store_true', help="use a single backbone for the context and feature encoders")
    parser.add_argument('--corr_levels', type=int, default=4, help="number of levels in the correlation pyramid")
    parser.add_argument('--corr_radius', type=int, default=4, help="width of the correlation pyramid")
    parser.add_argument('--n_downsample', type=int, default=2, help="resolution of the disparity field (1/2^K)")
    parser.add_argument('--context_norm', type=str, default="batch", choices=['group', 'batch', 'instance', 'none'], help="normalization of context encoder")
    parser.add_argument('--slow_fast_gru', action='store_true', help="iterate the low-res GRUs more frequently")
    parser.add_argument('--n_gru_layers', type=int, default=3, help="number of hidden GRU levels")
    parser.add_argument('-gpu','--ongpu',type=int, default=0, help="what gpu to use for inference")
    
    args = parser.parse_args()

    # print( args)
    print("Running demo with arguments:", args)

    demo(args)
