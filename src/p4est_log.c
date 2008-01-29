/* *INDENT-OFF* */

// Copyright (c) 2001, Bit Farm, Inc. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/* renamed from log4c.c to p4est_log.c and modified for p4est */

/**
 * Logging for C - implementation.
 *
 * See p4est_log.h for documentation.
 */

/*
 * the header p4est_log.h is included via p4est_base.h
 */
#include <p4est_base.h>

struct LogCategory _LOGV(LOG_ROOT_CAT) = {
    0, 0, 0,
    "root", P4EST_LP_UNINITIALIZED, 0,
    NULL, 0
};

static const char* s_controlString = NULL;

static const char* applyControlString(struct LogCategory* cat) {

    const char* cp = s_controlString;

    if (cp == NULL) return "Provide category with control string";

    while (*cp != 0) {
        const char *name, *dot, *eq;;
        cp += strspn(cp, " ");
        name = cp;
        cp += strcspn(cp, ".= ");
        dot = cp;
        cp += strcspn(cp, "= ");
        eq = cp;
        cp += strcspn(cp, " ");
        if (*dot != '.' || *eq != '=') {

            return "Invalid control string";
        }
        else if (0 == strncmp(name, cat->name, dot - name)) {
            if (0 == strncmp(dot + 1, "thresh", eq - dot - 1)) {
                log_setThreshold(cat, atoi(eq + 1));
            }
        }
    }
    return "";
} // applyControlString

void _log_logEvent(struct LogCategory* category, struct LogEvent* ev,
                   const char *fmt, ...)
{
    struct LogCategory* cat = category;
    va_start(ev->ap, fmt);
    while(1) {
        struct LogAppender* appender = cat->appender;
        if (appender != NULL) {
            appender->doAppend(appender, ev);
        }
        if (!cat->willLogToParent)
            break;

        cat = cat->parent;
    }
    va_end(ev->ap);
} // _log_logEvent

static const char * initCategory(struct LogCategory* category) {
    if (category == &_LOGV(LOG_ROOT_CAT)) {
        category->thresholdPriority = P4EST_LP_SILENT;
        category->appender = log_defaultLogAppender;
    } else {
        log_setParent(category, category->parent);
    }
    return applyControlString(category);
}

/**
 * This gets called the first time a category is referenced and performs the
 * initialization.
 * Also resets threshold to inherited!
 */
int _log_initCat(int priority, struct LogCategory* category) {

    initCategory(category);

    return priority >= category->thresholdPriority;
} // _log_initCat

void log_setParent(struct LogCategory* cat, struct LogCategory* parent) {

    P4EST_ASSERT (parent != NULL);

    // unlink from current parent
    if (cat->thresholdPriority != P4EST_LP_UNINITIALIZED) {
        struct LogCategory** cpp = &parent->firstChild;
        while(*cpp != cat && *cpp != NULL) {
            cpp = &(*cpp)->nextSibling;
        }
        P4EST_ASSERT (*cpp == cat);
        *cpp = cat->nextSibling;
    }

    // Set new parent
    cat->parent = parent;
    cat->nextSibling = parent->firstChild;
    parent->firstChild = cat;

    // Make sure parent is initialized
    if (parent->thresholdPriority == P4EST_LP_UNINITIALIZED) {
        initCategory(parent);
    }

    // Reset priority
    cat->thresholdPriority = parent->thresholdPriority;
    cat->isThreshInherited = 1;

} // log_setParent

static void setInheritedThresholds(struct LogCategory* cat)
{
    struct LogCategory* child = cat->firstChild;
    for( ; child != NULL; child = child->nextSibling)
    {
        if (child->isThreshInherited) {
            child->thresholdPriority = cat->thresholdPriority;
            setInheritedThresholds(child);
        }
    }
}

void log_setThreshold(struct LogCategory* cat, int thresholdPriority) {
    cat->thresholdPriority = thresholdPriority;
    cat->isThreshInherited = 0;
    setInheritedThresholds(cat);
}

const char* log_setControlString(const char* cs) {
    if (s_controlString == NULL) {
        s_controlString = cs;
        return initCategory(&_LOGV(LOG_ROOT_CAT));
    } else {
        return "log_setControlString should not be invoked twice.";
    }
}

void log_setAppender(struct LogCategory* cat, struct LogAppender* app) {
    cat->appender = app;
}

/**
 * The root category's default logging function.
 */

static char *priorityNames[] = {
    "Zero Priority",
    "TRACE",
    "DEBUG",
    "VERBOSE",
    "INFO",
    "STATISTICS",
    "PRODUCTION",
    "SILENT",
};

static void doAppend(struct LogAppender* this, struct LogEvent* ev);

static struct DefaultLogAppender {
    struct LogAppender appender;
    FILE *file;
} defaultLogAppender = { { doAppend }, NULL } ;

struct LogAppender* log_defaultLogAppender  = &defaultLogAppender.appender;

static void doAppend(struct LogAppender* this0, struct LogEvent* ev) {

    P4EST_CHECK_ABORT (0, "Default log appender should not be invoked");

    // TODO: define a format field in struct for timestamp, etc.
    char *pn;
    char buf[20];
    struct DefaultLogAppender* this = (struct DefaultLogAppender*)this0;
    const int spn = sizeof(priorityNames);

    if (this->file == NULL) this->file = stderr;

    if (ev->priority < 0) {
        pn = "Negative Priority NOT ALLOWED!!";
    }
    else if (ev->priority < spn) {
        pn = priorityNames[ev->priority];
    } else {
        sprintf(buf, "%s+%d",
                priorityNames[spn - 1], ev->priority - spn + 1);
    }
    fprintf(stderr, "%-7s ", pn);
    fprintf(stderr, "%s:%d: ", ev->fileName, ev->lineNum);
    vfprintf(stderr, ev->fmt, ev->ap);
    fprintf(stderr, "\n");
}

/* EOF p4est_log.c */
