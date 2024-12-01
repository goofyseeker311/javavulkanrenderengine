package fi.jkauppa.javavulkanclrenderengine;

import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.lwjgl.BufferUtils;
import org.lwjgl.PointerBuffer;
import org.lwjgl.glfw.Callbacks;
import org.lwjgl.glfw.GLFW;
import org.lwjgl.glfw.GLFWCursorPosCallbackI;
import org.lwjgl.glfw.GLFWErrorCallback;
import org.lwjgl.glfw.GLFWImage;
import org.lwjgl.glfw.GLFWImage.Buffer;
import org.lwjgl.glfw.GLFWKeyCallbackI;
import org.lwjgl.glfw.GLFWMouseButtonCallbackI;
import org.lwjgl.glfw.GLFWScrollCallbackI;
import org.lwjgl.glfw.GLFWVidMode;
import org.lwjgl.glfw.GLFWVulkan;
import org.lwjgl.openal.AL;
import org.lwjgl.openal.AL10;
import org.lwjgl.openal.ALC;
import org.lwjgl.openal.ALC11;
import org.lwjgl.openal.ALCCapabilities;
import org.lwjgl.openal.ALCapabilities;
import org.lwjgl.openal.EXTThreadLocalContext;
import org.lwjgl.system.Callback;
import org.lwjgl.system.MemoryStack;
import org.lwjgl.system.MemoryUtil;
import org.lwjgl.util.shaderc.Shaderc;
import org.lwjgl.util.shaderc.ShadercIncludeResolve;
import org.lwjgl.util.shaderc.ShadercIncludeResult;
import org.lwjgl.util.shaderc.ShadercIncludeResultRelease;
import org.lwjgl.vulkan.EXTDebugReport;
import org.lwjgl.vulkan.KHRDisplaySwapchain;
import org.lwjgl.vulkan.KHRSurface;
import org.lwjgl.vulkan.KHRSwapchain;
import org.lwjgl.vulkan.NVRayTracing;
import org.lwjgl.vulkan.VK13;
import org.lwjgl.vulkan.VkApplicationInfo;
import org.lwjgl.vulkan.VkAttachmentDescription;
import org.lwjgl.vulkan.VkAttachmentReference;
import org.lwjgl.vulkan.VkBufferCreateInfo;
import org.lwjgl.vulkan.VkCommandBuffer;
import org.lwjgl.vulkan.VkCommandBufferAllocateInfo;
import org.lwjgl.vulkan.VkCommandPoolCreateInfo;
import org.lwjgl.vulkan.VkDebugReportCallbackCreateInfoEXT;
import org.lwjgl.vulkan.VkDebugReportCallbackEXT;
import org.lwjgl.vulkan.VkDevice;
import org.lwjgl.vulkan.VkDeviceCreateInfo;
import org.lwjgl.vulkan.VkDeviceQueueCreateInfo;
import org.lwjgl.vulkan.VkGraphicsPipelineCreateInfo;
import org.lwjgl.vulkan.VkInstance;
import org.lwjgl.vulkan.VkInstanceCreateInfo;
import org.lwjgl.vulkan.VkLayerProperties;
import org.lwjgl.vulkan.VkMemoryAllocateInfo;
import org.lwjgl.vulkan.VkMemoryRequirements;
import org.lwjgl.vulkan.VkPhysicalDevice;
import org.lwjgl.vulkan.VkPhysicalDeviceMemoryProperties;
import org.lwjgl.vulkan.VkPipelineColorBlendAttachmentState;
import org.lwjgl.vulkan.VkPipelineColorBlendStateCreateInfo;
import org.lwjgl.vulkan.VkPipelineDepthStencilStateCreateInfo;
import org.lwjgl.vulkan.VkPipelineDynamicStateCreateInfo;
import org.lwjgl.vulkan.VkPipelineInputAssemblyStateCreateInfo;
import org.lwjgl.vulkan.VkPipelineLayoutCreateInfo;
import org.lwjgl.vulkan.VkPipelineMultisampleStateCreateInfo;
import org.lwjgl.vulkan.VkPipelineRasterizationStateCreateInfo;
import org.lwjgl.vulkan.VkPipelineShaderStageCreateInfo;
import org.lwjgl.vulkan.VkPipelineVertexInputStateCreateInfo;
import org.lwjgl.vulkan.VkPipelineViewportStateCreateInfo;
import org.lwjgl.vulkan.VkQueue;
import org.lwjgl.vulkan.VkQueueFamilyProperties;
import org.lwjgl.vulkan.VkRenderPassCreateInfo;
import org.lwjgl.vulkan.VkShaderModuleCreateInfo;
import org.lwjgl.vulkan.VkSubpassDependency;
import org.lwjgl.vulkan.VkSubpassDescription;
import org.lwjgl.vulkan.VkSurfaceFormatKHR;
import org.lwjgl.vulkan.VkVertexInputAttributeDescription;
import org.lwjgl.vulkan.VkVertexInputBindingDescription;

//import static org.lwjgl.vulkan.EXTDebugReport.*;
//import static org.lwjgl.vulkan.KHRSwapchain.*;
//import static org.lwjgl.vulkan.KHRSurface.*;
//import static org.lwjgl.vulkan.VK13.*;
//import static org.lwjgl.glfw.GLFWVulkan.*;

import static org.lwjgl.system.MemoryUtil.NULL;

import fi.jkauppa.javarenderengine.ModelLib;
import fi.jkauppa.javarenderengine.ModelLib.Entity;
import fi.jkauppa.javarenderengine.ModelLib.Sphere;
import fi.jkauppa.javarenderengine.ModelLib.Triangle;
import fi.jkauppa.javavulkanclrenderengine.ComputeLib;
import fi.jkauppa.javavulkanclrenderengine.ComputeLib.Device;
import fi.jkauppa.javarenderengine.UtilLib;

public class JavaVulkanCLRenderEngine {
	private static String programtitle = "Java OpenCL Render Engine v0.1.0.1";
	private int screenwidth = 0, screenheight = 0, graphicswidth = 0, graphicsheight = 0, graphicslength = 0;
	private boolean debug = System.getProperty("NDEBUG") == null;
    private String[] layers = {"VK_LAYER_LUNARG_standard_validation","VK_LAYER_KHRONOS_validation",};
	@SuppressWarnings("unused")
	private int litgraphicswidth = 0, litgraphicsheight = 0;
	private float graphicshfov = 70.0f, graphicsvfov = 39.375f;
	private long window = NULL;
	@SuppressWarnings("unused")
	private Callback debugProc = null;
	private int soundbuf = 0;
	private int sourcebuf = 0;
	private float frametime = 0.0f;
	private float frametimeavg = 0.0f;
	private ComputeLib computelib = null;
	private int selecteddevice = 0;
	@SuppressWarnings("unused")
	private boolean isfullscreen = false;
	private boolean vkinterop = true;
	private long opencldevice = NULL, openclqueue = NULL, openclprogram = NULL;
	private Device opencldevicedata = null;
	private String usingopencldevice = null;
	private long audiodevice = NULL;
	private long graphicsbufferptr = NULL, graphicszbufferptr = NULL, graphicshbufferptr = NULL, camposbufferptr = NULL, cammovbufferptr = NULL;
	private long tri1ptr = NULL, tri1lenptr = NULL, obj1ptr = NULL, obj1lenptr = NULL;
	private long tri2ptr = NULL, tri2lenptr = NULL, obj2ptr = NULL, obj2lenptr = NULL;
	private long trianglesptr = NULL, triangleslenptr = NULL, texturesptr = NULL, textureslenptr = NULL;
	@SuppressWarnings("unused")
	private long triangleslitptr = NULL;
	private long litptr = NULL;
	private long norptr = NULL;
	private float[] graphicsbuffer = null;
	@SuppressWarnings("unused")
	private float[] graphicszbuffer = null;
	private int[] graphicshbuffer = null;
	private float[] camerapos3fov2res2rotmat16 = null;
	private float[] cameramov3rot3 = null;
	private float[] trianglelist = null;
	private float[] trianglelist2 = null;
	private int[] trianglelistlength = {0};
	private int[] trianglelist2length = {0};
	private int[] triangleslistlen = {0};
	private int[] textureslist = null;
	private int[] textureslistlength = {0};
	private float[] objectlistpos3sca3rot3relsph4 = null;
	private float[] objectlist2pos3sca3rot3relsph4 = null;
	private int[] objectlistlength = {0};
	private int[] objectlist2length = {0};
	private int[] renderlit = {1};
	private int[] rendersphnorm = {0};
	private boolean keyfwd = false;
	private boolean keyback = false;
	private boolean keyleft = false;
	private boolean keyright = false;
	private boolean keyup = false;
	private boolean keydown = false;
	private boolean keyrleft = false;
	private boolean keyrright = false;
	private boolean keyspeed = false;
	private long nanolasttimetick = System.nanoTime();
	private double[] mousex = {0}, mousey = {0};
	private double lastmousex = 0, lastmousey = 0;
	private float lasttimedeltaseconds = 0.0f;
	private long monitor = NULL;
	private GLFWVidMode videomode = null;
	private KeyProcessor keyprocessor = new KeyProcessor();
	private MousePositionProcessor mouseposprocessor = new MousePositionProcessor();
	private MouseButtonProcessor mousebuttonprocessor = new MouseButtonProcessor();
	private MouseWheelProcessor mousewheelprocessor = new MouseWheelProcessor();

	public JavaVulkanCLRenderEngine(int vselecteddevice, int vfullscreen, int vvkinterop) {
		GLFWErrorCallback.createPrint(System.err).set();
		if (!GLFW.glfwInit()) {System.out.println("GLFW init failed."); System.exit(1);}
        if (!GLFWVulkan.glfwVulkanSupported()) {throw new AssertionError("GLFW vulkan init failed.");}
        PointerBuffer requiredExtensions = GLFWVulkan.glfwGetRequiredInstanceExtensions();
        if (requiredExtensions == null) {throw new AssertionError("GLFW vulkan extensions failed.");}
		this.monitor = GLFW.glfwGetPrimaryMonitor();
		this.videomode = GLFW.glfwGetVideoMode(this.monitor);
		this.screenwidth = 1280; this.screenheight = 720;
		long fullscreenmonitor = NULL;
		if (vfullscreen!=0) {
			this.isfullscreen = true;
			fullscreenmonitor = monitor;
			this.screenwidth = videomode.width();
			this.screenheight = videomode.height();
		}
		this.graphicswidth = screenwidth*2;
		this.graphicsheight = screenheight*2;
		if (vvkinterop==0) {
			this.vkinterop = false;
		}
		this.graphicshfov = (float)(Math.toDegrees(2.0f*Math.atan((((double)this.graphicswidth)/((double)this.graphicsheight))*Math.tan(Math.toRadians((double)(this.graphicsvfov/2.0f))))));
		this.graphicslength = this.graphicswidth*this.graphicsheight;
		GLFW.glfwDefaultWindowHints();
		//GLFW.glfwWindowHint(GLFW.GLFW_RED_BITS, videomode.redBits());
		//GLFW.glfwWindowHint(GLFW.GLFW_GREEN_BITS, videomode.greenBits());
		//GLFW.glfwWindowHint(GLFW.GLFW_BLUE_BITS, videomode.blueBits());
		GLFW.glfwWindowHint(GLFW.GLFW_CLIENT_API, GLFW.GLFW_NO_API);
		GLFW.glfwWindowHint(GLFW.GLFW_VISIBLE, GLFW.GLFW_FALSE);
		GLFW.glfwWindowHint(GLFW.GLFW_RESIZABLE, GLFW.GLFW_FALSE);
		if ((window=GLFW.glfwCreateWindow(screenwidth, screenheight, programtitle, fullscreenmonitor, NULL))==NULL) {System.out.println("GLFW create window failed."); System.exit(2);}
		GLFW.glfwSetInputMode(window, GLFW.GLFW_CURSOR, GLFW.GLFW_CURSOR_DISABLED);
		GLFW.glfwSetKeyCallback(window, keyprocessor);
		GLFW.glfwSetCursorPosCallback(window, mouseposprocessor);
		GLFW.glfwSetMouseButtonCallback(window, mousebuttonprocessor);
		GLFW.glfwSetScrollCallback(window, mousewheelprocessor);
		//GLFW.glfwMakeContextCurrent(window);
		//GLFW.glfwSwapInterval(1);
		GLFW.glfwShowWindow(window);
		//GLFW.glfwSwapBuffers(window);
		GLFW.glfwGetCursorPos(window, mousex, mousey);
		lastmousex = mousex[0]; lastmousey = mousey[0];
		BufferedImage iconimage = UtilLib.loadImage("res/images/icon.png", true);
		this.setIcon(iconimage);
		
		final VkInstance instance = createInstance(requiredExtensions);
        final VkDebugReportCallbackEXT debugCallback = new VkDebugReportCallbackEXT() {
            public int invoke(int flags, int objectType, long object, long location, int messageCode, long pLayerPrefix, long pMessage, long pUserData) {
                System.err.println("ERROR OCCURED: " + VkDebugReportCallbackEXT.getString(pMessage));
                return 0;
            }
        };
		final long debugCallbackHandle = setupDebugging(instance, EXTDebugReport.VK_DEBUG_REPORT_ERROR_BIT_EXT | EXTDebugReport.VK_DEBUG_REPORT_WARNING_BIT_EXT, debugCallback);
        final VkPhysicalDevice physicalDevice = getFirstPhysicalDevice(instance);
        final DeviceAndGraphicsQueueFamily deviceAndGraphicsQueueFamily = createDeviceAndGetGraphicsQueueFamily(physicalDevice);
        final VkDevice device = deviceAndGraphicsQueueFamily.device;
        int queueFamilyIndex = deviceAndGraphicsQueueFamily.queueFamilyIndex;
        final VkPhysicalDeviceMemoryProperties memoryProperties = deviceAndGraphicsQueueFamily.memoryProperties;
        LongBuffer pSurface = MemoryUtil.memAllocLong(1);
        int err = GLFWVulkan.glfwCreateWindowSurface(instance, window, null, pSurface);
        final long surface = pSurface.get(0);
        if (err != VK13.VK_SUCCESS) {throw new AssertionError("GLFW vulkan failed to create surface: " + translateVulkanResult(err));}
        
        final ColorFormatAndSpace colorFormatAndSpace = getColorFormatAndSpace(physicalDevice, surface);
        final long commandPool = createCommandPool(device, queueFamilyIndex);
        final VkCommandBuffer setupCommandBuffer = createCommandBuffer(device, commandPool);
        final VkQueue queue = createDeviceQueue(device, queueFamilyIndex);
        final long renderPass = createRenderPass(device, colorFormatAndSpace.colorFormat);
        final long renderCommandPool = createCommandPool(device, queueFamilyIndex);
        final Vertices vertices = createVertices(memoryProperties, device);
        final long pipeline = createPipeline(device, renderPass, vertices.createInfo);
        
		this.selecteddevice = vselecteddevice;
		this.computelib = new ComputeLib();
		this.opencldevice = this.computelib.devicelist[selecteddevice];
		this.opencldevicedata = this.computelib.devicemap.get(opencldevice);
		this.usingopencldevice = opencldevicedata.devicename;
		System.out.println("Using device["+selecteddevice+"]: "+opencldevicedata.devicename);
		this.openclqueue = opencldevicedata.queue;

		this.audiodevice = ALC11.alcOpenDevice((String)null);
		if (this.audiodevice == NULL) {
			throw new IllegalStateException("Failed to open default OpenAL device.");
		}
		ALCCapabilities audiodeviceCaps = ALC.createCapabilities(this.audiodevice);
		if (!audiodeviceCaps.OpenALC10) {
			throw new IllegalStateException();
		}
		long audiocontext = ALC11.alcCreateContext(this.audiodevice, (IntBuffer)null);
		boolean useTLC = audiodeviceCaps.ALC_EXT_thread_local_context && EXTThreadLocalContext.alcSetThreadContext(audiocontext);
		if (!useTLC) {
			if (!ALC11.alcMakeContextCurrent(audiocontext)) {
				throw new IllegalStateException();
			}
		}
		@SuppressWarnings("unused")
		ALCapabilities caps = AL.createCapabilities(audiodeviceCaps, MemoryUtil::memCallocPointer);
		this.soundbuf = AL10.alGenBuffers();
		this.sourcebuf = AL10.alGenSources();

		byte[] soundbytes = UtilLib.loadSound("res/sounds/firecannon.wav", 1, true);
		ByteBuffer soundbytesbuffer = MemoryUtil.memAlloc(soundbytes.length);
		soundbytesbuffer.put(soundbytes).rewind();
		AL10.alBufferData(this.soundbuf, AL10.AL_FORMAT_STEREO16, soundbytesbuffer, 44100);
		AL10.alSourcei(this.sourcebuf, AL10.AL_BUFFER, this.soundbuf);

		this.camerapos3fov2res2rotmat16 = new float[]{0.0f,0.0f,0.0f, graphicshfov,graphicsvfov, graphicswidth,graphicsheight, 1.0f,0.0f,0.0f,0.0f, 0.0f,1.0f,0.0f,0.0f, 0.0f,0.0f,1.0f,0.0f, 0.0f,0.0f,0.0f,1.0f};
		this.cameramov3rot3 = new float[]{0.0f,0.0f,0.0f, 0.0f,0.0f,0.0f};

		Entity loadmodel = ModelLib.loadOBJFileEntity("res/models/minef.obj", true);
		Entity loadmodel2 = ModelLib.loadOBJFileEntity("res/models/spaceboxgreen.obj", true);
		
		int imagecounter = 0;
		int imagelistlength = loadmodel.imagelist.length + loadmodel2.imagelist.length;
		
		BufferedImage textureimage = loadmodel2.imagelist[0];
		this.textureslistlength[0] = textureimage.getWidth();
		int texturesize = this.textureslistlength[0]*this.textureslistlength[0];
		this.textureslist = new int[imagelistlength*texturesize];
		for (int j=0;j<loadmodel.imagelist.length;j++) {
			textureimage = loadmodel.imagelist[j];
			DataBufferInt textureimagedataint = (DataBufferInt)textureimage.getRaster().getDataBuffer();
			int[] texturedata = textureimagedataint.getData();
			for (int i=0;i<texturesize;i++) {
				this.textureslist[(j+imagecounter)*texturesize+i] = texturedata[i];
			}
		}
		imagecounter += loadmodel.imagelist.length;
		for (int j=0;j<loadmodel2.imagelist.length;j++) {
			textureimage = loadmodel2.imagelist[j];
			DataBufferInt textureimagedataint = (DataBufferInt)textureimage.getRaster().getDataBuffer();
			int[] texturedata = textureimagedataint.getData();
			for (int i=0;i<texturesize;i++) {
				this.textureslist[(j+imagecounter)*texturesize+i] = texturedata[i];
			}
		}
		imagecounter += loadmodel2.imagelist.length;

		int imageoffset1 = 0;
		int imageoffset2 = loadmodel.imagelist.length;
		this.trianglelist = getEntityTriangles(loadmodel, imageoffset1);
		this.trianglelist2 = getEntityTriangles(loadmodel2, imageoffset2);
		this.trianglelistlength[0] = this.trianglelist.length/35;
		this.trianglelist2length[0] = this.trianglelist2.length/35;
		
		Sphere sphbv = loadmodel.sphereboundaryvolume;
		this.objectlistpos3sca3rot3relsph4 = new float[]{
				0.0f,0.0f,0.0f, 1.0f,1.0f,1.0f, 0.0f,0.0f,0.0f, (float)sphbv.x,(float)sphbv.y,(float)sphbv.z,(float)sphbv.r,
		};
		this.objectlistlength[0] = this.objectlistpos3sca3rot3relsph4.length/13;
		Sphere sphbv2 = loadmodel2.sphereboundaryvolume;
		this.objectlist2pos3sca3rot3relsph4 = new float[]{
				0.0f,0.0f,0.0f, 1.0f,1.0f,1.0f, 0.0f,0.0f,0.0f, (float)sphbv2.x,(float)sphbv2.y,(float)sphbv2.z,(float)sphbv2.r,
		};
		this.objectlist2length[0] = this.objectlist2pos3sca3rot3relsph4.length/13;
		
		this.triangleslistlen[0] = this.objectlistlength[0]*this.trianglelistlength[0] + this.objectlist2length[0]*this.trianglelist2length[0];
		this.litgraphicswidth = this.triangleslistlen[0] * 32;
		this.litgraphicsheight = 32*6;

		if (this.vkinterop) {
			//this.graphicsbufferptr = computelib.createSharedVKBuffer(opencldevice, buf);
		} else {
			this.graphicsbuffer = new float[graphicslength*4];
			this.graphicsbufferptr = computelib.createBuffer(opencldevice, graphicslength*4);
		}
		this.graphicszbufferptr = computelib.createBuffer(opencldevice, graphicslength);
		this.graphicszbuffer = new float[graphicslength];
		this.graphicshbufferptr = computelib.createBuffer(opencldevice, 1);
		this.graphicshbuffer = new int[1];
		this.camposbufferptr = computelib.createBuffer(opencldevice, camerapos3fov2res2rotmat16.length);
		computelib.writeBufferf(opencldevice, openclqueue, camposbufferptr, camerapos3fov2res2rotmat16);
		this.cammovbufferptr = computelib.createBuffer(opencldevice, cameramov3rot3.length);
		computelib.writeBufferf(opencldevice, openclqueue, cammovbufferptr, cameramov3rot3);

		this.trianglesptr = computelib.createBuffer(opencldevice, this.triangleslistlen[0]*35);
		this.triangleslenptr = computelib.createBuffer(opencldevice, 1);
		computelib.writeBufferi(opencldevice, openclqueue, triangleslenptr, this.triangleslistlen);
		this.triangleslitptr = computelib.createBuffer(opencldevice, this.triangleslistlen[0]*35);
		
		this.texturesptr = computelib.createBuffer(opencldevice, textureslist.length);
		computelib.writeBufferi(opencldevice, openclqueue, texturesptr, textureslist);
		this.textureslenptr = computelib.createBuffer(opencldevice, 1);
		computelib.writeBufferi(opencldevice, openclqueue, textureslenptr, textureslistlength);
		
		this.tri1ptr = computelib.createBuffer(opencldevice, trianglelist.length);
		computelib.writeBufferf(opencldevice, openclqueue, tri1ptr, trianglelist);
		this.tri1lenptr = computelib.createBuffer(opencldevice, 1);
		computelib.writeBufferi(opencldevice, openclqueue, tri1lenptr, trianglelistlength);
		this.obj1ptr = computelib.createBuffer(opencldevice, objectlistpos3sca3rot3relsph4.length);
		computelib.writeBufferf(opencldevice, openclqueue, obj1ptr, objectlistpos3sca3rot3relsph4);
		this.obj1lenptr = computelib.createBuffer(opencldevice, 1);
		computelib.writeBufferi(opencldevice, openclqueue, obj1lenptr, objectlistlength);

		this.tri2ptr = computelib.createBuffer(opencldevice, trianglelist2.length);
		computelib.writeBufferf(opencldevice, openclqueue, tri2ptr, trianglelist2);
		this.tri2lenptr = computelib.createBuffer(opencldevice, 1);
		computelib.writeBufferi(opencldevice, openclqueue, tri2lenptr, trianglelist2length);
		this.obj2ptr = computelib.createBuffer(opencldevice, objectlist2pos3sca3rot3relsph4.length);
		computelib.writeBufferf(opencldevice, openclqueue, obj2ptr, objectlist2pos3sca3rot3relsph4);
		this.obj2lenptr = computelib.createBuffer(opencldevice, 1);
		computelib.writeBufferi(opencldevice, openclqueue, obj2lenptr, objectlist2length);
		
		this.litptr = computelib.createBuffer(opencldevice, 1);
		computelib.writeBufferi(opencldevice, openclqueue, litptr, renderlit);
		this.norptr = computelib.createBuffer(opencldevice, 1);
		computelib.writeBufferi(opencldevice, openclqueue, norptr, rendersphnorm);
		
		String programSource = new String(ComputeLib.loadProgram("res/clprograms/programlib.cl", true));
		this.openclprogram = this.computelib.compileProgram(opencldevice, programSource);
		System.out.println("init.");
	}

	public void run() {
		while(!GLFW.glfwWindowShouldClose(window)) {
			long nanonewtimetick = System.nanoTime();
			lasttimedeltaseconds = (nanonewtimetick - nanolasttimetick)/1000000000.0f;
			nanolasttimetick = nanonewtimetick;
			tick(lasttimedeltaseconds);
			if (this.vkinterop) {computelib.acquireSharedBuffer(openclqueue, graphicsbufferptr);}
			render();
			if (this.vkinterop) {computelib.releaseSharedBuffer(openclqueue, graphicsbufferptr);}
			if (!this.vkinterop) {}
			//GLFW.glfwSwapBuffers(window);
			GLFW.glfwPollEvents();
		}

		Callbacks.glfwFreeCallbacks(window);
		GLFW.glfwDestroyWindow(window);
		GLFW.glfwTerminate();
		GLFW.glfwSetErrorCallback(null).free();
	}

	public static void main(String[] args) {
		System.out.println(programtitle);
		int argdevice = 0;
		int argfullscreen = 1;
		int argglinterop = 1;
		try {argdevice = Integer.parseInt(args[0]);} catch(Exception ex) {}
		try {argfullscreen = Integer.parseInt(args[1]);} catch(Exception ex) {}
		try {argglinterop = Integer.parseInt(args[2]);} catch(Exception ex) {}
		JavaVulkanCLRenderEngine app = new JavaVulkanCLRenderEngine(argdevice, argfullscreen, argglinterop);
		app.run();
		System.out.println("exit.");
		System.exit(0);
	}

	private void tick(float deltatimeseconds) {
		float ds = deltatimeseconds;
		GLFW.glfwSetWindowTitle(window, programtitle+": "+String.format("%.0f",1000.0f/frametimeavg).replace(',', '.')+
				"fps, computetime: "+String.format("%.3f",frametimeavg).replace(',', '.')+"ms ["+usingopencldevice+"] ("
				+screenwidth+"x"+screenheight+") tickdeltatime: "+String.format("%.0f",deltatimeseconds*1000.0f)+"ms"
				+" ["+(this.vkinterop?"VKINTEROP":"COPYBUFFER")+"]"
				);
		cameramov3rot3[0] = 0.0f;
		cameramov3rot3[1] = 0.0f;
		cameramov3rot3[2] = 0.0f;
		cameramov3rot3[3] = 0.0f;
		cameramov3rot3[4] = 0.0f;
		cameramov3rot3[5] = 0.0f;
		float sp = this.keyspeed?100.0f:1.0f;
		if (this.keyfwd) {cameramov3rot3[0] = ds*sp;}
		if (this.keyback) {cameramov3rot3[0] = -ds*sp;}
		if (this.keyleft) {cameramov3rot3[1] = -ds*sp;}
		if (this.keyright) {cameramov3rot3[1] = ds*sp;}
		if (this.keyup) {cameramov3rot3[2] = ds*sp;}
		if (this.keydown) {cameramov3rot3[2] = -ds*sp;}
		if (this.keyrleft) {cameramov3rot3[5] = -ds;}
		if (this.keyrright) {cameramov3rot3[5] = ds;}
		cameramov3rot3[4] = -(float)(0.001f*(mousex[0]-lastmousex));
		cameramov3rot3[3] = (float)(0.001f*(mousey[0]-lastmousey));
		lastmousex = mousex[0];
		lastmousey = mousey[0];
	}

	public void render() {
		long framestarttime = System.nanoTime();
		computelib.writeBufferf(opencldevice, openclqueue, cammovbufferptr, cameramov3rot3);
		computelib.runProgram(opencldevice, openclqueue, openclprogram, "movecamerakernel", new long[]{camposbufferptr,cammovbufferptr}, new int[]{0}, new int[]{1});
		computelib.insertBarrier(openclqueue);
		computelib.readBufferf(opencldevice, openclqueue, camposbufferptr, camerapos3fov2res2rotmat16);
		objectlist2pos3sca3rot3relsph4[0] = camerapos3fov2res2rotmat16[0];
		objectlist2pos3sca3rot3relsph4[1] = camerapos3fov2res2rotmat16[1];
		objectlist2pos3sca3rot3relsph4[2] = camerapos3fov2res2rotmat16[2];
		computelib.writeBufferf(opencldevice, openclqueue, obj2ptr, objectlist2pos3sca3rot3relsph4);
		computelib.runProgram(opencldevice, openclqueue, openclprogram, "clearviewkernel", new long[]{graphicsbufferptr,graphicszbufferptr,graphicshbufferptr,camposbufferptr}, new int[]{0,0}, new int[]{graphicswidth,2});
		int trianglecount1 = this.objectlistlength[0]*this.trianglelistlength[0]*35;
		computelib.runProgram(opencldevice, openclqueue, openclprogram, "transformobjectkernel", new long[]{trianglesptr,tri1ptr,tri1lenptr,obj1ptr,obj1lenptr}, new int[]{0,0,0}, new int[]{1,objectlistlength[0],trianglelistlength[0]});
		computelib.runProgram(opencldevice, openclqueue, openclprogram, "transformobjectkernel", new long[]{trianglesptr,tri2ptr,tri2lenptr,obj2ptr,obj2lenptr}, new int[]{trianglecount1,0,0}, new int[]{1,objectlist2length[0],trianglelist2length[0]});
		computelib.insertBarrier(openclqueue);
		//computelib.runProgram(opencldevice, openclqueue, program, "lightobjectkernel", new long[]{,,,triangleslitptr,trianglesptr,triangleslenptr,texturesptr,textureslenptr}, new int[]{0,0,0}, new int[]{triangleslistlen[0],1,triangleslistlen[0]});
		computelib.insertBarrier(openclqueue);
		computelib.runProgram(opencldevice, openclqueue, openclprogram, "renderplaneviewkernel", new long[]{graphicsbufferptr,graphicszbufferptr,graphicshbufferptr,camposbufferptr,trianglesptr,triangleslenptr,texturesptr,textureslenptr,litptr,norptr}, new int[]{0,0}, new int[]{graphicswidth,2});
		computelib.insertBarrier(openclqueue);
		computelib.runProgram(opencldevice, openclqueue, openclprogram, "rendercrosskernel", new long[]{graphicsbufferptr,graphicszbufferptr,graphicshbufferptr,camposbufferptr}, new int[]{0}, new int[]{1});
		computelib.waitForQueue(openclqueue);
		computelib.readBufferi(opencldevice, openclqueue, graphicshbufferptr, graphicshbuffer);
		if (!this.vkinterop) {
			float[] newgraphicsbuffer = new float[graphicslength*4];
			computelib.readBufferf(opencldevice, openclqueue, graphicsbufferptr, newgraphicsbuffer);
			graphicsbuffer = newgraphicsbuffer;
		}
		long frameendtime = System.nanoTime();
		frametime = (frameendtime-framestarttime)/1000000.0f;
		frametimeavg = frametimeavg*0.9f+frametime*0.1f;
	}

    private VkInstance createInstance(PointerBuffer requiredExtensions) {
        VkApplicationInfo appInfo = VkApplicationInfo.calloc()
                .sType$Default()
                .apiVersion(VK13.VK_API_VERSION_1_3);
        PointerBuffer ppEnabledExtensionNames = MemoryUtil.memAllocPointer(requiredExtensions.remaining() + 1);
        ppEnabledExtensionNames.put(requiredExtensions);
        ByteBuffer VK_EXT_DEBUG_REPORT_EXTENSION = MemoryUtil.memUTF8(EXTDebugReport.VK_EXT_DEBUG_REPORT_EXTENSION_NAME);
        ppEnabledExtensionNames.put(VK_EXT_DEBUG_REPORT_EXTENSION);
        ppEnabledExtensionNames.flip();
        PointerBuffer ppEnabledLayerNames = debug ? allocateLayerBuffer(layers) : null;
        VkInstanceCreateInfo pCreateInfo = VkInstanceCreateInfo.calloc().sType$Default().pApplicationInfo(appInfo).ppEnabledExtensionNames(ppEnabledExtensionNames).ppEnabledLayerNames(ppEnabledLayerNames);
        PointerBuffer pInstance = MemoryUtil.memAllocPointer(1);
        int err = VK13.vkCreateInstance(pCreateInfo, null, pInstance);
        long instance = pInstance.get(0);
        MemoryUtil.memFree(pInstance);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to create VkInstance: " + translateVulkanResult(err));
        }
        VkInstance ret = new VkInstance(instance, pCreateInfo);
        pCreateInfo.free();
        if(ppEnabledLayerNames != null) MemoryUtil.memFree(ppEnabledLayerNames);
        MemoryUtil.memFree(VK_EXT_DEBUG_REPORT_EXTENSION);
        MemoryUtil.memFree(ppEnabledExtensionNames);
        MemoryUtil.memFree(appInfo.pApplicationName());
        MemoryUtil.memFree(appInfo.pEngineName());
        appInfo.free();
        return ret;
    }

    public String translateVulkanResult(int result) {
        switch (result) {
        case VK13.VK_SUCCESS:
            return "Command successfully completed.";
        case VK13.VK_NOT_READY:
            return "A fence or query has not yet completed.";
        case VK13.VK_TIMEOUT:
            return "A wait operation has not completed in the specified time.";
        case VK13.VK_EVENT_SET:
            return "An event is signaled.";
        case VK13.VK_EVENT_RESET:
            return "An event is unsignaled.";
        case VK13.VK_INCOMPLETE:
            return "A return array was too small for the result.";
        case KHRSwapchain.VK_SUBOPTIMAL_KHR:
            return "A swapchain no longer matches the surface properties exactly, but can still be used to present to the surface successfully.";

        case VK13.VK_ERROR_OUT_OF_HOST_MEMORY:
            return "A host memory allocation has failed.";
        case VK13.VK_ERROR_OUT_OF_DEVICE_MEMORY:
            return "A device memory allocation has failed.";
        case VK13.VK_ERROR_INITIALIZATION_FAILED:
            return "Initialization of an object could not be completed for implementation-specific reasons.";
        case VK13.VK_ERROR_DEVICE_LOST:
            return "The logical or physical device has been lost.";
        case VK13.VK_ERROR_MEMORY_MAP_FAILED:
            return "Mapping of a memory object has failed.";
        case VK13.VK_ERROR_LAYER_NOT_PRESENT:
            return "A requested layer is not present or could not be loaded.";
        case VK13.VK_ERROR_EXTENSION_NOT_PRESENT:
            return "A requested extension is not supported.";
        case VK13.VK_ERROR_FEATURE_NOT_PRESENT:
            return "A requested feature is not supported.";
        case VK13.VK_ERROR_INCOMPATIBLE_DRIVER:
            return "The requested version of Vulkan is not supported by the driver or is otherwise incompatible for implementation-specific reasons.";
        case VK13.VK_ERROR_TOO_MANY_OBJECTS:
            return "Too many objects of the type have already been created.";
        case VK13.VK_ERROR_FORMAT_NOT_SUPPORTED:
            return "A requested format is not supported on this device.";
        case KHRSurface.VK_ERROR_SURFACE_LOST_KHR:
            return "A surface is no longer available.";
        case KHRSurface.VK_ERROR_NATIVE_WINDOW_IN_USE_KHR:
            return "The requested window is already connected to a VkSurfaceKHR, or to some other non-Vulkan API.";
        case KHRSwapchain.VK_ERROR_OUT_OF_DATE_KHR:
            return "A surface has changed in such a way that it is no longer compatible with the swapchain, and further presentation requests using the "
                    + "swapchain will fail. Applications must query the new surface properties and recreate their swapchain if they wish to continue"
                    + "presenting to the surface.";
        case KHRDisplaySwapchain.VK_ERROR_INCOMPATIBLE_DISPLAY_KHR:
            return "The display used by a swapchain does not use the same presentable image layout, or is incompatible in a way that prevents sharing an"
                    + " image.";
        case EXTDebugReport.VK_ERROR_VALIDATION_FAILED_EXT:
            return "A validation layer found an error.";
        default:
            return String.format("%s [%d]", "Unknown", Integer.valueOf(result));
        }
    }

    private long setupDebugging(VkInstance instance, int flags, VkDebugReportCallbackEXT callback) {
        VkDebugReportCallbackCreateInfoEXT dbgCreateInfo = VkDebugReportCallbackCreateInfoEXT.calloc()
                .sType$Default()
                .pfnCallback(callback)
                .flags(flags);
        LongBuffer pCallback = MemoryUtil.memAllocLong(1);
        int err = EXTDebugReport.vkCreateDebugReportCallbackEXT(instance, dbgCreateInfo, null, pCallback);
        long callbackHandle = pCallback.get(0);
        MemoryUtil.memFree(pCallback);
        dbgCreateInfo.free();
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to create VkInstance: " + translateVulkanResult(err));
        }
        return callbackHandle;
    }
    
    public final PointerBuffer allocateLayerBuffer(String[] layers) {
        final Set<String> availableLayers = getAvailableLayers();
        PointerBuffer ppEnabledLayerNames = MemoryUtil.memAllocPointer(layers.length);
        for (int i = 0; i < layers.length; i++) {
            final String layer = layers[i];
            if (availableLayers.contains(layer)) {
                ppEnabledLayerNames.put(MemoryUtil.memUTF8(layer));
            }
        }
        ppEnabledLayerNames.flip();
        return ppEnabledLayerNames;
    }

    private final Set<String> getAvailableLayers() {
        final Set<String> res = new HashSet<>();
        final int[] ip = new int[1];
        VK13.vkEnumerateInstanceLayerProperties(ip, null);
        final int count = ip[0];

        try (final MemoryStack stack = MemoryStack.stackPush()) {
            if (count > 0) {
                final VkLayerProperties.Buffer instanceLayers = VkLayerProperties.malloc(count, stack);
                VK13.vkEnumerateInstanceLayerProperties(ip, instanceLayers);
                for (int i = 0; i < count; i++) {
                    final String layerName = instanceLayers.get(i).layerNameString();
                    res.add(layerName);
                }
            }
        }

        return res;
    }
    
    private class DeviceAndGraphicsQueueFamily {
        VkDevice device;
        int queueFamilyIndex;
        VkPhysicalDeviceMemoryProperties memoryProperties;
    }

    private VkPhysicalDevice getFirstPhysicalDevice(VkInstance instance) {
        IntBuffer pPhysicalDeviceCount = MemoryUtil.memAllocInt(1);
        int err = VK13.vkEnumeratePhysicalDevices(instance, pPhysicalDeviceCount, null);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to get number of physical devices: " + translateVulkanResult(err));
        }
        PointerBuffer pPhysicalDevices = MemoryUtil.memAllocPointer(pPhysicalDeviceCount.get(0));
        err = VK13.vkEnumeratePhysicalDevices(instance, pPhysicalDeviceCount, pPhysicalDevices);
        long physicalDevice = pPhysicalDevices.get(0);
        MemoryUtil.memFree(pPhysicalDeviceCount);
        MemoryUtil.memFree(pPhysicalDevices);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to get physical devices: " + translateVulkanResult(err));
        }
        return new VkPhysicalDevice(physicalDevice, instance);
    }

    private DeviceAndGraphicsQueueFamily createDeviceAndGetGraphicsQueueFamily(VkPhysicalDevice physicalDevice) {
        IntBuffer pQueueFamilyPropertyCount = MemoryUtil.memAllocInt(1);
        VK13.vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, pQueueFamilyPropertyCount, null);
        int queueCount = pQueueFamilyPropertyCount.get(0);
        VkQueueFamilyProperties.Buffer queueProps = VkQueueFamilyProperties.calloc(queueCount);
        VK13.vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, pQueueFamilyPropertyCount, queueProps);
        MemoryUtil.memFree(pQueueFamilyPropertyCount);
        int graphicsQueueFamilyIndex;
        for (graphicsQueueFamilyIndex = 0; graphicsQueueFamilyIndex < queueCount; graphicsQueueFamilyIndex++) {
            if ((queueProps.get(graphicsQueueFamilyIndex).queueFlags() & VK13.VK_QUEUE_GRAPHICS_BIT) != 0)
                break;
        }
        queueProps.free();
        FloatBuffer pQueuePriorities = MemoryUtil.memAllocFloat(1).put(0.0f);
        pQueuePriorities.flip();
        VkDeviceQueueCreateInfo.Buffer queueCreateInfo = VkDeviceQueueCreateInfo.calloc(1).sType$Default().queueFamilyIndex(graphicsQueueFamilyIndex).pQueuePriorities(pQueuePriorities);

        PointerBuffer extensions = MemoryUtil.memAllocPointer(1);
        ByteBuffer VK_KHR_SWAPCHAIN_EXTENSION = MemoryUtil.memUTF8(KHRSwapchain.VK_KHR_SWAPCHAIN_EXTENSION_NAME);
        extensions.put(VK_KHR_SWAPCHAIN_EXTENSION);
        extensions.flip();

        VkDeviceCreateInfo deviceCreateInfo = VkDeviceCreateInfo.calloc().sType$Default().pQueueCreateInfos(queueCreateInfo).ppEnabledExtensionNames(extensions);

        PointerBuffer pDevice = MemoryUtil.memAllocPointer(1);
        int err = VK13.vkCreateDevice(physicalDevice, deviceCreateInfo, null, pDevice);
        long device = pDevice.get(0);
        MemoryUtil.memFree(pDevice);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to create device: " + translateVulkanResult(err));
        }

        VkPhysicalDeviceMemoryProperties memoryProperties = VkPhysicalDeviceMemoryProperties.calloc();
        VK13.vkGetPhysicalDeviceMemoryProperties(physicalDevice, memoryProperties);

        DeviceAndGraphicsQueueFamily ret = new DeviceAndGraphicsQueueFamily();
        ret.device = new VkDevice(device, physicalDevice, deviceCreateInfo);
        ret.queueFamilyIndex = graphicsQueueFamilyIndex;
        ret.memoryProperties = memoryProperties;

        deviceCreateInfo.free();
        MemoryUtil.memFree(VK_KHR_SWAPCHAIN_EXTENSION);
        MemoryUtil.memFree(extensions);
        MemoryUtil.memFree(pQueuePriorities);
        return ret;
    }

    private class ColorFormatAndSpace {
        int colorFormat;
        int colorSpace;
    }
    
    private ColorFormatAndSpace getColorFormatAndSpace(VkPhysicalDevice physicalDevice, long surface) {
        IntBuffer pQueueFamilyPropertyCount = MemoryUtil.memAllocInt(1);
        VK13.vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, pQueueFamilyPropertyCount, null);
        int queueCount = pQueueFamilyPropertyCount.get(0);
        VkQueueFamilyProperties.Buffer queueProps = VkQueueFamilyProperties.calloc(queueCount);
        VK13.vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, pQueueFamilyPropertyCount, queueProps);
        MemoryUtil.memFree(pQueueFamilyPropertyCount);
        IntBuffer supportsPresent = MemoryUtil.memAllocInt(queueCount);
        for (int i = 0; i < queueCount; i++) {
            supportsPresent.position(i);
            int err = KHRSurface.vkGetPhysicalDeviceSurfaceSupportKHR(physicalDevice, i, surface, supportsPresent);
            if (err != VK13.VK_SUCCESS) {
                throw new AssertionError("Failed to physical device surface support: " + translateVulkanResult(err));
            }
        }
        int graphicsQueueNodeIndex = Integer.MAX_VALUE;
        int presentQueueNodeIndex = Integer.MAX_VALUE;
        for (int i = 0; i < queueCount; i++) {
            if ((queueProps.get(i).queueFlags() & VK13.VK_QUEUE_GRAPHICS_BIT) != 0) {
                if (graphicsQueueNodeIndex == Integer.MAX_VALUE) {
                    graphicsQueueNodeIndex = i;
                }
                if (supportsPresent.get(i) == VK13.VK_TRUE) {
                    graphicsQueueNodeIndex = i;
                    presentQueueNodeIndex = i;
                    break;
                }
            }
        }
        queueProps.free();
        if (presentQueueNodeIndex == Integer.MAX_VALUE) {
            for (int i = 0; i < queueCount; ++i) {
                if (supportsPresent.get(i) == VK13.VK_TRUE) {
                    presentQueueNodeIndex = i;
                    break;
                }
            }
        }
        MemoryUtil.memFree(supportsPresent);
        if (graphicsQueueNodeIndex == Integer.MAX_VALUE) {throw new AssertionError("No graphics queue found");}
        if (presentQueueNodeIndex == Integer.MAX_VALUE) {throw new AssertionError("No presentation queue found");}
        if (graphicsQueueNodeIndex != presentQueueNodeIndex) {throw new AssertionError("Presentation queue != graphics queue");}
        IntBuffer pFormatCount = MemoryUtil.memAllocInt(1);
        int err = KHRSurface.vkGetPhysicalDeviceSurfaceFormatsKHR(physicalDevice, surface, pFormatCount, null);
        int formatCount = pFormatCount.get(0);
        if (err != VK13.VK_SUCCESS) {throw new AssertionError("Failed to query number of physical device surface formats: " + translateVulkanResult(err));}

        VkSurfaceFormatKHR.Buffer surfFormats = VkSurfaceFormatKHR.calloc(formatCount);
        err = KHRSurface.vkGetPhysicalDeviceSurfaceFormatsKHR(physicalDevice, surface, pFormatCount, surfFormats);
        MemoryUtil.memFree(pFormatCount);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to query physical device surface formats: " + translateVulkanResult(err));
        }

        int colorFormat;
        if (formatCount == 1 && surfFormats.get(0).format() == VK13.VK_FORMAT_UNDEFINED) {
            colorFormat = VK13.VK_FORMAT_B8G8R8A8_UNORM;
        } else {
            colorFormat = surfFormats.get(0).format();
        }
        int colorSpace = surfFormats.get(0).colorSpace();
        surfFormats.free();

        ColorFormatAndSpace ret = new ColorFormatAndSpace();
        ret.colorFormat = colorFormat;
        ret.colorSpace = colorSpace;
        return ret;
    }

    private long createCommandPool(VkDevice device, int queueNodeIndex) {
        VkCommandPoolCreateInfo cmdPoolInfo = VkCommandPoolCreateInfo.calloc()
                .sType$Default()
                .queueFamilyIndex(queueNodeIndex)
                .flags(VK13.VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT);
        LongBuffer pCmdPool = MemoryUtil.memAllocLong(1);
        int err = VK13.vkCreateCommandPool(device, cmdPoolInfo, null, pCmdPool);
        long commandPool = pCmdPool.get(0);
        cmdPoolInfo.free();
        MemoryUtil.memFree(pCmdPool);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to create command pool: " + translateVulkanResult(err));
        }
        return commandPool;
    }

    private VkQueue createDeviceQueue(VkDevice device, int queueFamilyIndex) {
        PointerBuffer pQueue = MemoryUtil.memAllocPointer(1);
        VK13.vkGetDeviceQueue(device, queueFamilyIndex, 0, pQueue);
        long queue = pQueue.get(0);
        MemoryUtil.memFree(pQueue);
        return new VkQueue(queue, device);
    }

    private VkCommandBuffer createCommandBuffer(VkDevice device, long commandPool) {
        VkCommandBufferAllocateInfo cmdBufAllocateInfo = VkCommandBufferAllocateInfo.calloc()
                .sType$Default()
                .commandPool(commandPool)
                .level(VK13.VK_COMMAND_BUFFER_LEVEL_PRIMARY)
                .commandBufferCount(1);
        PointerBuffer pCommandBuffer = MemoryUtil.memAllocPointer(1);
        int err = VK13.vkAllocateCommandBuffers(device, cmdBufAllocateInfo, pCommandBuffer);
        cmdBufAllocateInfo.free();
        long commandBuffer = pCommandBuffer.get(0);
        MemoryUtil.memFree(pCommandBuffer);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to allocate command buffer: " + translateVulkanResult(err));
        }
        return new VkCommandBuffer(commandBuffer, device);
    }

    private long createRenderPass(VkDevice device, int colorFormat) {
        VkAttachmentDescription.Buffer attachments = VkAttachmentDescription.calloc(1).format(colorFormat).samples(VK13.VK_SAMPLE_COUNT_1_BIT).loadOp(VK13.VK_ATTACHMENT_LOAD_OP_CLEAR)
                .storeOp(VK13.VK_ATTACHMENT_STORE_OP_STORE).stencilLoadOp(VK13.VK_ATTACHMENT_LOAD_OP_DONT_CARE).stencilStoreOp(VK13.VK_ATTACHMENT_STORE_OP_DONT_CARE)
                .initialLayout(VK13.VK_IMAGE_LAYOUT_UNDEFINED).finalLayout(KHRSwapchain.VK_IMAGE_LAYOUT_PRESENT_SRC_KHR);
        VkAttachmentReference.Buffer colorReference = VkAttachmentReference.calloc(1).attachment(0).layout(VK13.VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL);
        VkSubpassDescription.Buffer subpass = VkSubpassDescription.calloc(1).pipelineBindPoint(VK13.VK_PIPELINE_BIND_POINT_GRAPHICS)
                .colorAttachmentCount(colorReference.remaining()).pColorAttachments(colorReference);
        VkSubpassDependency.Buffer dependency = VkSubpassDependency.calloc(1).srcSubpass(VK13.VK_SUBPASS_EXTERNAL).srcStageMask(VK13.VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT)
                .dstAccessMask(VK13.VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT).dstStageMask(VK13.VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT).dependencyFlags(VK13.VK_DEPENDENCY_BY_REGION_BIT);
        VkRenderPassCreateInfo renderPassInfo = VkRenderPassCreateInfo.calloc().sType$Default().pAttachments(attachments).pSubpasses(subpass).pDependencies(dependency);

        LongBuffer pRenderPass = MemoryUtil.memAllocLong(1);
        int err = VK13.vkCreateRenderPass(device, renderPassInfo, null, pRenderPass);
        long renderPass = pRenderPass.get(0);
        MemoryUtil.memFree(pRenderPass);
        dependency.free();
        renderPassInfo.free();
        colorReference.free();
        subpass.free();
        attachments.free();
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to create clear render pass: " + translateVulkanResult(err));
        }
        return renderPass;
    }

    private boolean getMemoryType(VkPhysicalDeviceMemoryProperties deviceMemoryProperties, int typeBits, int properties, IntBuffer typeIndex) {
        int bits = typeBits;
        for (int i = 0; i < 32; i++) {
            if ((bits & 1) == 1) {
                if ((deviceMemoryProperties.memoryTypes(i).propertyFlags() & properties) == properties) {
                    typeIndex.put(0, i);
                    return true;
                }
            }
            bits >>= 1;
        }
        return false;
    }
    
    private class Vertices {
        long verticesBuf;
        VkPipelineVertexInputStateCreateInfo createInfo;
    }

    private Vertices createVertices(VkPhysicalDeviceMemoryProperties deviceMemoryProperties, VkDevice device) {
        ByteBuffer vertexBuffer = MemoryUtil.memAlloc(3 * 2 * 4);
        FloatBuffer fb = vertexBuffer.asFloatBuffer();
        fb.put(-0.5f).put(-0.5f);
        fb.put( 0.5f).put(-0.5f);
        fb.put( 0.0f).put( 0.5f);

        VkMemoryAllocateInfo memAlloc = VkMemoryAllocateInfo.calloc().sType$Default();
        VkMemoryRequirements memReqs = VkMemoryRequirements.calloc();

        int err;
        VkBufferCreateInfo bufInfo = VkBufferCreateInfo.calloc()
                .sType$Default()
                .size(vertexBuffer.remaining())
                .usage(VK13.VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
        LongBuffer pBuffer = MemoryUtil.memAllocLong(1);
        err = VK13.vkCreateBuffer(device, bufInfo, null, pBuffer);
        long verticesBuf = pBuffer.get(0);
        MemoryUtil.memFree(pBuffer);
        bufInfo.free();
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to create vertex buffer: " + translateVulkanResult(err));
        }

        VK13.vkGetBufferMemoryRequirements(device, verticesBuf, memReqs);
        memAlloc.allocationSize(memReqs.size());
        IntBuffer memoryTypeIndex = MemoryUtil.memAllocInt(1);
        getMemoryType(deviceMemoryProperties, memReqs.memoryTypeBits(), VK13.VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT, memoryTypeIndex);
        memAlloc.memoryTypeIndex(memoryTypeIndex.get(0));
        MemoryUtil.memFree(memoryTypeIndex);
        memReqs.free();

        LongBuffer pMemory = MemoryUtil.memAllocLong(1);
        err = VK13.vkAllocateMemory(device, memAlloc, null, pMemory);
        long verticesMem = pMemory.get(0);
        MemoryUtil.memFree(pMemory);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to allocate vertex memory: " + translateVulkanResult(err));
        }

        PointerBuffer pData = MemoryUtil.memAllocPointer(1);
        err = VK13.vkMapMemory(device, verticesMem, 0, memAlloc.allocationSize(), 0, pData);
        memAlloc.free();
        long data = pData.get(0);
        MemoryUtil.memFree(pData);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to map vertex memory: " + translateVulkanResult(err));
        }

        MemoryUtil.memCopy(MemoryUtil.memAddress(vertexBuffer), data, vertexBuffer.remaining());
        MemoryUtil.memFree(vertexBuffer);
        VK13.vkUnmapMemory(device, verticesMem);
        err = VK13.vkBindBufferMemory(device, verticesBuf, verticesMem, 0);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to bind memory to vertex buffer: " + translateVulkanResult(err));
        }

        VkVertexInputBindingDescription.Buffer bindingDescriptor = VkVertexInputBindingDescription.calloc(1).binding(0).stride(2 * 4).inputRate(VK13.VK_VERTEX_INPUT_RATE_VERTEX);
        VkVertexInputAttributeDescription.Buffer attributeDescriptions = VkVertexInputAttributeDescription.calloc(1);
        attributeDescriptions.get(0).binding(0).location(0).format(VK13.VK_FORMAT_R32G32_SFLOAT).offset(0);
        VkPipelineVertexInputStateCreateInfo vi = VkPipelineVertexInputStateCreateInfo.calloc().sType$Default().pVertexBindingDescriptions(bindingDescriptor).pVertexAttributeDescriptions(attributeDescriptions);

        Vertices ret = new Vertices();
        ret.createInfo = vi;
        ret.verticesBuf = verticesBuf;
        return ret;
    }

    private long createPipeline(VkDevice device, long renderPass, VkPipelineVertexInputStateCreateInfo vi) {
        int err;
        VkPipelineInputAssemblyStateCreateInfo inputAssemblyState = VkPipelineInputAssemblyStateCreateInfo.calloc().sType$Default().topology(VK13.VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST);
        VkPipelineRasterizationStateCreateInfo rasterizationState = VkPipelineRasterizationStateCreateInfo.calloc().sType$Default().polygonMode(VK13.VK_POLYGON_MODE_FILL)
                .cullMode(VK13.VK_CULL_MODE_NONE).frontFace(VK13.VK_FRONT_FACE_COUNTER_CLOCKWISE).lineWidth(1.0f);
        VkPipelineColorBlendAttachmentState.Buffer colorWriteMask = VkPipelineColorBlendAttachmentState.calloc(1).colorWriteMask(0xF);
        VkPipelineColorBlendStateCreateInfo colorBlendState = VkPipelineColorBlendStateCreateInfo.calloc().sType$Default().pAttachments(colorWriteMask);
        VkPipelineViewportStateCreateInfo viewportState = VkPipelineViewportStateCreateInfo.calloc().sType$Default().viewportCount(1).scissorCount(1);

        IntBuffer pDynamicStates = MemoryUtil.memAllocInt(2);
        pDynamicStates.put(VK13.VK_DYNAMIC_STATE_VIEWPORT).put(VK13.VK_DYNAMIC_STATE_SCISSOR).flip();
        VkPipelineDynamicStateCreateInfo dynamicState = VkPipelineDynamicStateCreateInfo.calloc().sType$Default().pDynamicStates(pDynamicStates);
        VkPipelineDepthStencilStateCreateInfo depthStencilState = VkPipelineDepthStencilStateCreateInfo.calloc().sType$Default().depthCompareOp(VK13.VK_COMPARE_OP_ALWAYS);
        depthStencilState.back().failOp(VK13.VK_STENCIL_OP_KEEP).passOp(VK13.VK_STENCIL_OP_KEEP).compareOp(VK13.VK_COMPARE_OP_ALWAYS);
        depthStencilState.front(depthStencilState.back());

        VkPipelineMultisampleStateCreateInfo multisampleState = VkPipelineMultisampleStateCreateInfo.calloc().sType$Default().rasterizationSamples(VK13.VK_SAMPLE_COUNT_1_BIT);
        VkPipelineShaderStageCreateInfo.Buffer shaderStages = VkPipelineShaderStageCreateInfo.calloc(2);
        shaderStages.get(0).set(loadShader(device, "res/vkshaders/triangle.vert", VK13.VK_SHADER_STAGE_VERTEX_BIT));
        shaderStages.get(1).set(loadShader(device, "res/vkshaders/triangle.frag", VK13.VK_SHADER_STAGE_FRAGMENT_BIT));

        VkPipelineLayoutCreateInfo pPipelineLayoutCreateInfo = VkPipelineLayoutCreateInfo.calloc().sType$Default();

        LongBuffer pPipelineLayout = MemoryUtil.memAllocLong(1);
        err = VK13.vkCreatePipelineLayout(device, pPipelineLayoutCreateInfo, null, pPipelineLayout);
        long layout = pPipelineLayout.get(0);
        MemoryUtil.memFree(pPipelineLayout);
        pPipelineLayoutCreateInfo.free();
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to create pipeline layout: " + translateVulkanResult(err));
        }

        VkGraphicsPipelineCreateInfo.Buffer pipelineCreateInfo = VkGraphicsPipelineCreateInfo.calloc(1).sType$Default().layout(layout).renderPass(renderPass)
                .pVertexInputState(vi).pInputAssemblyState(inputAssemblyState).pRasterizationState(rasterizationState).pColorBlendState(colorBlendState)
                .pMultisampleState(multisampleState).pViewportState(viewportState).pDepthStencilState(depthStencilState).pStages(shaderStages).pDynamicState(dynamicState);

        LongBuffer pPipelines = MemoryUtil.memAllocLong(1);
        err = VK13.vkCreateGraphicsPipelines(device, VK13.VK_NULL_HANDLE, pipelineCreateInfo, null, pPipelines);
        long pipeline = pPipelines.get(0);
        shaderStages.free();
        multisampleState.free();
        depthStencilState.free();
        dynamicState.free();
        MemoryUtil.memFree(pDynamicStates);
        viewportState.free();
        colorBlendState.free();
        colorWriteMask.free();
        rasterizationState.free();
        inputAssemblyState.free();
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to create pipeline: " + translateVulkanResult(err));
        }
        return pipeline;
    }

    private long loadShader(String classPath, VkDevice device, int stage) {
        ByteBuffer shaderCode = glslToSpirv(classPath, stage);
        int err;
        VkShaderModuleCreateInfo moduleCreateInfo = VkShaderModuleCreateInfo.calloc()
                .sType$Default()
                .pCode(shaderCode);
        LongBuffer pShaderModule = MemoryUtil.memAllocLong(1);
        err = VK13.vkCreateShaderModule(device, moduleCreateInfo, null, pShaderModule);
        long shaderModule = pShaderModule.get(0);
        MemoryUtil.memFree(pShaderModule);
        if (err != VK13.VK_SUCCESS) {
            throw new AssertionError("Failed to create shader module: " + translateVulkanResult(err));
        }
        return shaderModule;
    }
    
    private VkPipelineShaderStageCreateInfo loadShader(VkDevice device, String classPath, int stage) {
        VkPipelineShaderStageCreateInfo shaderStage = VkPipelineShaderStageCreateInfo.calloc()
                .sType$Default()
                .stage(stage)
                .module(loadShader(classPath, device, stage))
                .pName(MemoryUtil.memUTF8("main"));
        return shaderStage;
    }
    
    public ByteBuffer glslToSpirv(String classPath, int vulkanStage) {
    	byte[] sourceShader = ComputeLib.loadProgram(classPath, true);
		ByteBuffer src = BufferUtils.createByteBuffer(sourceShader.length);
		src.put(sourceShader).rewind();
        long compiler = Shaderc.shaderc_compiler_initialize();
        long options = Shaderc.shaderc_compile_options_initialize();
        ShadercIncludeResolve resolver;
        ShadercIncludeResultRelease releaser;
        Shaderc.shaderc_compile_options_set_target_env(options, Shaderc.shaderc_target_env_vulkan, Shaderc.shaderc_env_version_vulkan_1_3);
        Shaderc.shaderc_compile_options_set_target_spirv(options, Shaderc.shaderc_spirv_version_1_4);
        Shaderc.shaderc_compile_options_set_optimization_level(options, Shaderc.shaderc_optimization_level_performance);
        Shaderc.shaderc_compile_options_set_include_callbacks(options, resolver = new ShadercIncludeResolve() {
            public long invoke(long user_data, long requested_source, int type, long requesting_source, long include_depth) {
                ShadercIncludeResult res = ShadercIncludeResult.calloc();
                String src = classPath.substring(0, classPath.lastIndexOf('/')) + "/" + MemoryUtil.memUTF8(requested_source);
                byte[] srcShader = ComputeLib.loadProgram(src, true);
                if (srcShader!=null) {
	        		ByteBuffer srcbytes = BufferUtils.createByteBuffer(srcShader.length);
	        		srcbytes.put(srcShader).rewind();
	                res.content(srcbytes);
	                res.source_name(MemoryUtil.memUTF8(src));
	                return res.address();
                } else {
                    throw new AssertionError("Failed to resolve include: " + src);
                }
            }
        }, releaser = new ShadercIncludeResultRelease() {
            public void invoke(long user_data, long include_result) {
                ShadercIncludeResult result = ShadercIncludeResult.create(include_result);
                MemoryUtil.memFree(result.source_name());
                result.free();
            }
        }, 0L);
        long res;
        try (MemoryStack stack = MemoryStack.stackPush()) {
            res = Shaderc.shaderc_compile_into_spv(compiler, src, vulkanStageToShadercKind(vulkanStage), stack.UTF8(classPath), stack.UTF8("main"), options);
            if (res == 0L)
                throw new AssertionError("Internal error during compilation!");
        }
        if (Shaderc.shaderc_result_get_compilation_status(res) != Shaderc.shaderc_compilation_status_success) {
            throw new AssertionError("Shader compilation failed: " + Shaderc.shaderc_result_get_error_message(res));
        }
        int size = (int) Shaderc.shaderc_result_get_length(res);
        ByteBuffer resultBytes = BufferUtils.createByteBuffer(size);
        resultBytes.put(Shaderc.shaderc_result_get_bytes(res));
        resultBytes.flip();
        Shaderc.shaderc_result_release(res);
        Shaderc.shaderc_compiler_release(compiler);
        releaser.free();
        resolver.free();
        return resultBytes;
    }
    
    private int vulkanStageToShadercKind(int stage) {
        switch (stage) {
        case VK13.VK_SHADER_STAGE_VERTEX_BIT:
            return Shaderc.shaderc_vertex_shader;
        case VK13.VK_SHADER_STAGE_FRAGMENT_BIT:
            return Shaderc.shaderc_fragment_shader;
        case NVRayTracing.VK_SHADER_STAGE_RAYGEN_BIT_NV:
            return Shaderc.shaderc_raygen_shader;
        case NVRayTracing.VK_SHADER_STAGE_CLOSEST_HIT_BIT_NV:
            return Shaderc.shaderc_closesthit_shader;
        case NVRayTracing.VK_SHADER_STAGE_MISS_BIT_NV:
            return Shaderc.shaderc_miss_shader;
        case NVRayTracing.VK_SHADER_STAGE_ANY_HIT_BIT_NV:
            return Shaderc.shaderc_anyhit_shader;
        case NVRayTracing.VK_SHADER_STAGE_INTERSECTION_BIT_NV:
            return Shaderc.shaderc_intersection_shader;
        case VK13.VK_SHADER_STAGE_COMPUTE_BIT:
            return Shaderc.shaderc_compute_shader;
        default:
            throw new IllegalArgumentException("Stage: " + stage);
        }
    }
    
	private void setIcon(BufferedImage iconimage) {
		DataBufferInt iconimagedataint = (DataBufferInt)iconimage.getRaster().getDataBuffer();
		int[] iconimageints = iconimagedataint.getData();
		IntBuffer iconimageintbuffer = IntBuffer.wrap(iconimageints);
		ByteBuffer iconimagebytebuffer = MemoryUtil.memAlloc(iconimageints.length*4);
		iconimagebytebuffer.asIntBuffer().put(iconimageintbuffer);
		for (int i=0;i<iconimageints.length;i++) {
			byte cr = iconimagebytebuffer.get(i*4+2);
			byte cg = iconimagebytebuffer.get(i*4+1);
			byte cb = iconimagebytebuffer.get(i*4+0);
			byte ca = iconimagebytebuffer.get(i*4+3);
			iconimagebytebuffer.put(i*4+0, cr);
			iconimagebytebuffer.put(i*4+1, cg);
			iconimagebytebuffer.put(i*4+2, cb);
			iconimagebytebuffer.put(i*4+3, ca);
		}
		Buffer iconimagebuffer = GLFWImage.create(1);
		GLFWImage iconglfwimage = GLFWImage.create().set(iconimage.getWidth(), iconimage.getHeight(), iconimagebytebuffer);
		iconimagebuffer.put(0, iconglfwimage);
		GLFW.glfwSetWindowIcon(window, iconimagebuffer);
		MemoryUtil.memFree(iconimagebytebuffer);
	}
	
	private float[] getEntityTriangles(Entity loadmodel, int imageidoffset) {
		ArrayList<Float> trianglearraylist = new ArrayList<Float>();
		for (int j=0;j<loadmodel.childlist.length;j++) {
			Entity object = loadmodel.childlist[j];
			for (int i=0;i<object.trianglelist.length;i++) {
				Triangle modeltri = object.trianglelist[i];
				trianglearraylist.add((float)modeltri.pos1.x);
				trianglearraylist.add((float)modeltri.pos1.y);
				trianglearraylist.add((float)modeltri.pos1.z);
				trianglearraylist.add((float)modeltri.pos2.x);
				trianglearraylist.add((float)modeltri.pos2.y);
				trianglearraylist.add((float)modeltri.pos2.z);
				trianglearraylist.add((float)modeltri.pos3.x);
				trianglearraylist.add((float)modeltri.pos3.y);
				trianglearraylist.add((float)modeltri.pos3.z);
				trianglearraylist.add((float)modeltri.norm.dx);
				trianglearraylist.add((float)modeltri.norm.dy);
				trianglearraylist.add((float)modeltri.norm.dz);
				trianglearraylist.add((float)modeltri.pos1.tex.u);
				trianglearraylist.add((float)modeltri.pos1.tex.v);
				trianglearraylist.add((float)modeltri.pos2.tex.u);
				trianglearraylist.add((float)modeltri.pos2.tex.v);
				trianglearraylist.add((float)modeltri.pos3.tex.u);
				trianglearraylist.add((float)modeltri.pos3.tex.v);
				trianglearraylist.add((float)modeltri.mat.imageid+imageidoffset);
				float[] matfacecolor = modeltri.mat.facecolor.getRGBComponents(new float[4]);
				trianglearraylist.add((float)matfacecolor[0]);
				trianglearraylist.add((float)matfacecolor[1]);
				trianglearraylist.add((float)matfacecolor[2]);
				trianglearraylist.add((float)matfacecolor[3]);
				float[] matemissivecolor = modeltri.mat.emissivecolor.getRGBComponents(new float[4]);
				trianglearraylist.add((float)matemissivecolor[0]);
				trianglearraylist.add((float)matemissivecolor[1]);
				trianglearraylist.add((float)matemissivecolor[2]);
				trianglearraylist.add((float)matemissivecolor[3]);
				float[] lightmapcolor = {0.0f,0.0f,0.0f,0.0f};
				if (modeltri.mat.ambientcolor!=null) {lightmapcolor = modeltri.mat.ambientcolor.getRGBComponents(new float[4]);}
				trianglearraylist.add((float)lightmapcolor[0]);
				trianglearraylist.add((float)lightmapcolor[1]);
				trianglearraylist.add((float)lightmapcolor[2]);
				trianglearraylist.add((float)lightmapcolor[3]);
				trianglearraylist.add((float)modeltri.mat.roughness);
				trianglearraylist.add((float)modeltri.mat.metallic);
				trianglearraylist.add((float)modeltri.mat.refraction);
				trianglearraylist.add((float)modeltri.mat.transparency);
			}
		}

		Float[] trianglefloats= trianglearraylist.toArray(new Float[trianglearraylist.size()]);
		float[] trianglelist = new float[trianglefloats.length];
		for (int i=0;i<trianglefloats.length;i++) {
			trianglelist[i] = trianglefloats[i];
		}
		return trianglelist;
	}

	private class KeyProcessor implements GLFWKeyCallbackI {
		@Override public void invoke(long window, int key, int scancode, int action, int mods) {
			if (action==GLFW.GLFW_PRESS) {
				if (key==GLFW.GLFW_KEY_W) {keyfwd = true;}
				if (key==GLFW.GLFW_KEY_S) {keyback = true;}
				if (key==GLFW.GLFW_KEY_A) {keyleft = true;}
				if (key==GLFW.GLFW_KEY_D) {keyright = true;}
				if (key==GLFW.GLFW_KEY_SPACE) {keyup = true;}
				if (key==GLFW.GLFW_KEY_LEFT_SHIFT) {keydown = true;}
				if (key==GLFW.GLFW_KEY_Q) {keyrleft = true;}
				if (key==GLFW.GLFW_KEY_E) {keyrright = true;}
				if (key==GLFW.GLFW_KEY_LEFT_CONTROL) {keyspeed = true;}
			}
			if (action==GLFW.GLFW_RELEASE) {
				if (key==GLFW.GLFW_KEY_W) {keyfwd = false;}
				if (key==GLFW.GLFW_KEY_S) {keyback = false;}
				if (key==GLFW.GLFW_KEY_A) {keyleft = false;}
				if (key==GLFW.GLFW_KEY_D) {keyright = false;}
				if (key==GLFW.GLFW_KEY_SPACE) {keyup = false;}
				if (key==GLFW.GLFW_KEY_LEFT_SHIFT) {keydown = false;}
				if (key==GLFW.GLFW_KEY_Q) {keyrleft = false;}
				if (key==GLFW.GLFW_KEY_E) {keyrright = false;}
				if (key==GLFW.GLFW_KEY_LEFT_CONTROL) {keyspeed = false;}
			}
		}
	}
	private class MousePositionProcessor implements GLFWCursorPosCallbackI {
		@Override public void invoke(long window, double xpos, double ypos) {
			mousex[0] = xpos;
			mousey[0] = ypos;
		}
	}
	private class MouseButtonProcessor implements GLFWMouseButtonCallbackI {
		@Override public void invoke(long window, int button, int action, int mods) {
			if ((button==0)&&(action==1)) {
				AL10.alSourcePlay(sourcebuf);
			}
		}
	}
	private class MouseWheelProcessor implements GLFWScrollCallbackI {
		@Override public void invoke(long window, double xoffset, double yoffset) {
			System.out.println("xoffset: "+xoffset+" yoffset: "+yoffset);
		}
	}


}
