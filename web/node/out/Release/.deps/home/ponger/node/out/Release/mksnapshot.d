cmd_/home/ponger/node/out/Release/mksnapshot := flock /home/ponger/node/out/Release/linker.lock g++ -pthread -rdynamic -m64 -m64  -o /home/ponger/node/out/Release/mksnapshot -Wl,--start-group /home/ponger/node/out/Release/obj.target/mksnapshot/deps/v8/src/mksnapshot.o /home/ponger/node/out/Release/obj.target/deps/v8/tools/gyp/libv8_base.a /home/ponger/node/out/Release/obj.target/deps/v8/tools/gyp/libv8_nosnapshot.a /home/ponger/node/out/Release/obj.target/deps/v8/tools/gyp/libv8_libplatform.a /home/ponger/node/out/Release/obj.target/deps/v8/tools/gyp/libv8_libbase.a -Wl,--end-group -lrt